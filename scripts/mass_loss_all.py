#!/usr/bin/env python
# coding: utf-8

# In[1]:


import xarray as xr
import matplotlib.pyplot as plt
import os
import sys
import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.spatial import KDTree
import cmocean.cm as ocm


proj_dir = os.path.join(os.pardir)
data_dir = os.path.join(proj_dir,'data')

src_dir = os.path.join(proj_dir,'src')
sys.path.append(src_dir)

from tools.log_progress import log_progress
from visualization.shiftedColorMap import shiftedColorMap

coords_path = os.path.join(data_dir,'analysis','external','antbounds','shelves2.mat')
#coords_path = os.path.join(data_dir,'analysis','external','antbounds','shelves10.mat')
int_dir = os.path.join(data_dir,'analysis','interim')
pro_dir = os.path.join(data_dir,'analysis','processed')
fig_dir = os.path.join(os.pardir,os.pardir,'reports','figures')
grid_path = os.path.join(data_dir,'preprocessing','processed','waom2_grd.nc')
#grid_path = os.path.join(data_dir,'preprocessing','processed','waom10_grd.nc')
#m_path = os.path.join(data_dir,'analysis','raw','waom4','ocean_avg_0009.nc')
#temp_path = os.path.join(data_dir,'analysis','raw','waom2','ocean_avg_0538-0610_temp_avg.nc')
#zeta_path = os.path.join(data_dir,'analysis','raw','waom2','ocean_avg_0538-0610_zeta_avg.nc')
#avg_path = os.path.join(data_dir,'analysis','raw','waom10','ocean_avg_0006.nc')

pd.options.display.float_format = '{:,.2f}'.format


# In[16]:


from dask.distributed import Client
c = Client()
c


# In[3]:


grd = xr.open_dataset(grid_path)
m_path = os.path.join(data_dir,'analysis','raw','waom2','ocean_avg_0538-0610_m_avg.nc')
m = xr.open_dataset(m_path).m.isel(ocean_time=0)

#m = xr.open_dataset(m_path).m.squeeze()
#temp = xr.open_mfdataset(temp_path).temp.squeeze()
#zeta = xr.open_dataset(zeta_path).zeta.squeeze()
#m = avg.m
#temp = avg.temp
#zeta = avg.zeta

s2a = 3600*24*365.25
rhoi = 916
mask_vostock = (grd.lat_rho<-75) & (grd.lat_rho>-80) & (grd.lon_rho>95) & (grd.lon_rho<115)


# # Creating masks
# ## Create ice shelf masks using ant_bounds coordinates

# In[5]:


# find nearest neighbours of ant_bounds lat_lon coords on the roms grid and 
# define ice shelf masks and middle coordinates
shelves = {}

# read in ice shelf coordinates derived from ant_bounds matlab package
coords = sio.loadmat(coords_path)['shelves'][0].squeeze()

for idx in range(np.size(coords)):
    name = coords[idx][0][0]
    shelves[name]={}
    shelves[name]['lat']= coords[idx][1].squeeze()
    shelves[name]['lon']= coords[idx][2].squeeze() 

lat_flat = grd.lat_rho.stack(etaxi = ('eta_rho','xi_rho'))
lon_flat = grd.lon_rho.stack(etaxi = ('eta_rho','xi_rho'))
points = np.column_stack((lat_flat.values,lon_flat.values))
tree = KDTree(points)

for name,data in log_progress(shelves.items(),name='ice shelves'):
    
    lats = data['lat']
    lons = data['lon']

    target = np.column_stack((lats,lons))
    dist, ind = tree.query(target,k=1)
    
    eta = lat_flat[np.unique(ind)].eta_rho.values
    xi = lat_flat[np.unique(ind)].xi_rho.values
    shelves[name]['eta'] = eta
    shelves[name]['xi'] = xi
    
    mask_tmp = np.zeros_like(grd.mask_rho.values)
    mask_tmp[eta,xi] = 1
    mask_tmp[grd.zice.values == 0.0] = 0
    mask_tmp[grd.mask_rho.values == 0.0] = 0
    mask_tmp = (mask_tmp == 1)  
    shelves[name]['mask']=mask_tmp
    
    xi_m,eta_m = np.median(xi),np.median(eta)
    shelves[name]['xi_mid']=xi_m
    shelves[name]['eta_mid']=eta_m


# In[8]:


shelves = {k:v for k,v in shelves.items() if v['mask'].any()}
dict_out_path = os.path.join(int_dir,'shelves2_masks.npy')
np.save(dict_out_path,shelves)


# ## create size masks

# In[4]:


shelves_path = os.path.join(data_dir,'analysis','interim','shelves2_masks.npy')
shelves = np.load(shelves_path,allow_pickle=True).item()


# In[5]:


dA = (grd.pm*grd.pn)**-1/10**9
for k,v in shelves.items():
    shelves[k]['A'] = dA.where(v['mask']).sum()


# In[6]:


def make_mask_size(low_lim,up_lim):
    
    shelves_sel = {name:data for (name,data) in shelves.items() if ((data['A'] >= low_lim) & (data['A'] < up_lim))}

    mask = np.zeros_like(grd.mask_rho,dtype=bool)

    for k,v in shelves_sel.items():
        mask[v['mask']] = True
        
    return xr.DataArray(mask,dims=('eta_rho','xi_rho'))

mask_size = {}
mask_size['all'] = (grd.mask_rho == 1 ) & (grd.zice < 0.0)
mask_size['small'] = make_mask_size(0,5)
mask_size['medium'] = make_mask_size(5,45)
mask_size['large'] = make_mask_size(45,10000)
mask_size['tiny'] = mask_size['all'] & ~mask_size['small'] & ~mask_size['medium'] & ~mask_size['large'] & ~(mask_vostock)


# In[153]:


dict_out_path = os.path.join(int_dir,'mask_size2.npy')
np.save(dict_out_path,mask_size)


# In[7]:


size=np.zeros_like(grd.mask_rho.values)

for i,mask_name in zip([1,2,3,4],['tiny','small','medium','large']):
    size[mask_size[mask_name]] = i
size = xr.DataArray(size,dims=('eta_rho','xi_rho'))   
size = size.where((grd.mask_rho ==1) & (grd.zice<0.0))


# In[23]:


mask_ice = (grd.mask_rho ==1) & (grd.zice < 0) & (mask_vostock == 0)
plt.close()
fig,ax = plt.subplots(figsize = (15,10))
ax.contourf(grd.mask_rho.where((grd.mask_rho==0) | (mask_vostock == 1)),colors=(('0.6','0.6','0.6')))
#colors=['blue','orange','green','red']
img = size.plot.contourf(ax=ax,levels=[0,1,2,3,4],cmap='viridis',add_colorbar=False)
ax.contourf(grd.mask_rho.where((mask_vostock == 1)),colors=(('0.6','0.6','0.6')))
#grd.mask_rho.where((mask_ice == 1) & (wb_rho == 0)).plot.contourf(ax=ax, colors=('gray'),add_colorbar=False)
for key,data in log_progress(shelves.items(),name='shelves'):
           ax.contour(data['mask'], colors=('black'),linewidths=0.05)
ax.set_aspect('equal')
ax.axis('off')
cbar = plt.colorbar(img, ticks=[0.5,1.5,2.5,3.5])
cbar.ax.tick_params(labelsize=18)
cbar.set_label('Area ($10^3$ km$^2$)',size=18)
cbar.ax.set_yticklabels(['not in \nMEaSURES','< 5','5 - 45',' > 45'])
plt.tight_layout()
plt.savefig(os.path.join(fig_dir,'size_mask.png'),transparent=True,dpi=300)
plt.show()


# ## create mask sector

# In[10]:


mask_sector = {}
mask_sector['Total Antarctica'] = xr.DataArray(np.ones_like(grd.lon_rho,dtype=bool),dims=('eta_rho','xi_rho'))
mask_sector['Western East Antarctica'] = (grd.lon_rho>=-10.0) & (grd.lon_rho<60)
mask_sector['Amery/Prydz Bay'] = (grd.lon_rho>=60.0) & (grd.lon_rho<80.0)
mask_sector['Sabrina Coast/Aurora subglacial basin'] = (grd.lon_rho>=80.0) & (grd.lon_rho<130.0) & ~(mask_vostock)
mask_sector['George V Coast/Wilkes subglacial basin'] = ((grd.lon_rho>=130.0) & (grd.lon_rho<155.0)) | ((grd.lon_rho>=155.0) & (grd.lon_rho<170.0) & (grd.lat_rho>=-72))
mask_sector['Ross Sea'] = (grd.lon_rho>=155.0) & ~mask_sector['George V Coast/Wilkes subglacial basin'] |(grd.lon_rho<-140.0)|((grd.lon_rho>=-140.0) & (grd.lon_rho<-120.0) &(grd.lat_rho<-77.0)) 
mask_sector['Amundsen Sea'] = ((grd.lon_rho>=-140.0) & (grd.lon_rho<-120.0) & (grd.lat_rho>=-77.0)) |((grd.lon_rho>=-120.0) & (grd.lon_rho<-90.0))
mask_sector['Bellingshausen Sea'] = (grd.lon_rho>=-90.0) & (grd.lon_rho<-66.0) & (grd.lat_rho>=-75.0)
mask_sector['Weddell Sea'] = (grd.lon_rho > -90) & (grd.lon_rho < -10) & ~ mask_sector['Bellingshausen Sea']


# In[15]:


dict_out_path = os.path.join(int_dir,'mask_sector2.npy')
np.save(dict_out_path,mask_sector)


# In[66]:


plt.close()
fig,ax = plt.subplots()
for k,v in mask_sector.items():
    ax.contour(v,colors=('black'),linewidths=0.01)
ax.set_xticks([])
ax.set_yticks([])
ax.set_aspect('equal')
plt.savefig(os.path.join(fig_dir,'tmp','lines.png'),transparent=True,dpi=300)
plt.show()


# ## create mask sector on the continental shelf

# In[4]:


mask_shelf = ((grd.h < 1000) & (grd.mask_rho ==1) & (grd.zice == 0.0)) | ((grd.zice < 0.0) & (grd.mask_rho ==1))
mask_shelf_sector = {}
for key,item in mask_sector.items():
    mask_shelf_sector[key] = (item & mask_shelf)
mask_shelf_sector['Amundsen Sea'] &= grd.lat_rho < -70 
mask_shelf_sector['Ross Sea'] &= grd.lat_rho < -70
mask_shelf_sector['George V Coast/Wilkes subglacial basin'] &= ~((grd.lon_rho > 160) & (grd.lat_rho > -68.5)) 
mask_shelf_sector['Sabrina Coast/Aurora subglacial basin'] &= grd.lat_rho <= -64
mask_shelf_sector['Western East Antarctica'] &= grd.lat_rho < -64
mask_shelf_sector['Weddell Sea'] &= ~((grd.lat_rho > -63) & (grd.lon_rho > -48))

mask_shelf_sector['Total Antarctica'] = mask_shelf_sector['Amundsen Sea'].copy()
for key,item in mask_shelf_sector.items():
    mask_shelf_sector['Total Antarctica'] |= item


# In[6]:


plt.close()
mask_shelf_sector['Total Antarctica'].plot()
plt.show()


# In[7]:


dict_out_path = os.path.join(int_dir,'mask_shelf_sector2.npy')
np.save(dict_out_path,mask_shelf_sector)


# ## create mask front

# In[67]:


def make_mask_front(grd,nb_cells):
    
    mask_rho = grd.mask_rho.values
    mask_land = np.zeros_like(mask_rho)
    mask_land[mask_rho == 0] = 1
    mask_zice = np.zeros_like(mask_land)
    mask_zice[grd.zice.values*mask_rho != 0] = 1

    mask_front = np.zeros_like(grd.mask_rho.values)

    for j in grd.eta_rho.values:
        for i in grd. xi_rho.values:
            if mask_zice[j,i] == 1:
                j_min = max(j-nb_cells,0)
                j_max = min(j+nb_cells, np.size(mask_rho,0))
                i_min = max(i-nb_cells,0)
                i_max = min(i+nb_cells+1, np.size(mask_rho,1))

                if np.any(mask_zice[j_min:j_max,i_min:i_max] + mask_land[j_min:j_max,i_min:i_max]== 0):
                        mask_front[j,i] = 1
                        
    grd['mask_front'] = (('eta_rho','xi_rho'),mask_front)
    
    return grd
grd = make_mask_front(grd,3)


# In[56]:


get_ipython().run_line_magic('matplotlib', 'notebook')
levels = [-200]
mask_shelf = (grd.mask_rho == 1) & (grd.zice < 0.0)
plt.close()
fig,ax = plt.subplots(figsize=(8,8))
(avg.where(mask_shelf).m*s2a).plot(ax=ax,vmin=-2,vmax=2,cmap='bwr')
cs = grd.zice.where(grd.mask_rho).plot.contour(levels=levels,ax=ax,colors='k')
grd.mask_front.where(grd.mask_rho).plot.contour(levels=[0.1],ax=ax,colors='k',linestyles='dashed')
ax.set_aspect('equal')
plt.show()


# ## depth mask

# In[11]:


get_ipython().run_line_magic('matplotlib', 'inline')

depth1 = -200
depth2 = -400


depth_mask = grd.mask_rho.copy()
depth_mask.values[grd.zice == 0.0] = 0

depth_mask.values[(grd.zice >= depth1) & (grd.zice < 0.0)] = 1
depth_mask.values[(grd.zice >= depth2) & (grd.zice < depth1)] = 2
depth_mask.values[grd.zice < depth2] = 3

plt.close()
fig,ax = plt.subplots(figsize=(12,8))
img = depth_mask.where(grd.mask_rho).plot(ax=ax,add_colorbar=False)
plt.colorbar(img, ax=ax)
ax.set_aspect('equal')
plt.show()


# In[18]:


mask_ice = (grd.mask_rho == 1) & (grd.zice < 0.0) & (mask_vostock== 0)
mask_depth = {}
mask_depth['all'] = mask_ice
mask_depth['<200'] = (mask_ice == 1) & (grd.zice>=-200)
mask_depth['<400'] =  (mask_ice == 1) & (grd.zice>=-400)
mask_depth['>400'] = (mask_ice == 1) & (grd.zice<-400)


# In[10]:

mask_sector_path = os.path.join(int_dir,'mask_sector2.npy') 
mask_sector = np.load(mask_sector_path,allow_pickle=True).item()

mask_depth_path = os.path.join(int_dir,'mask_depth_2.npy')
mask_depth = np.load(mask_depth_path).item()


# # Analysis
# ## mass loss for different depth ranges

# In[28]:


mask = (grd.mask_rho==1) & (grd.zice < 0.0)
y = m.values[mask]*s2a
x = grd.zice.values[mask]*-1
plt.close()
fig,ax = plt.subplots()
plt.plot(x,y,'k.')
plt.xlabel('depth in m')
plt.ylabel('ismr in m/yr')
plt.savefig(os.path.join(fig_dir,'ismr_vs_depth.png'),transparent=True,dpi=300)
plt.show()


# In[29]:


bins = np.arange(0,2800,100)
depths = grd.zice*-1
dA = (1/(grd.pm*grd.pn)).where(mask)
ismr2bmb = dA*rhoi*(10**-12)

ismr = (m.where(mask).groupby_bins(depths,bins).mean()*s2a).to_series()
A = dA.groupby_bins(depths,bins).sum().to_series()/10**9
bmb =  ((m.where(mask)*ismr2bmb).groupby_bins(depths,bins).sum()*s2a).to_series()


# In[76]:


plt.rcParams.update({'font.size': 18})
plt.close()
fig,axes = plt.subplots(3,figsize=(12,12))
ax1,ax2,ax3 = axes.flatten()
ismr.plot.bar(ax=ax1,color='k')
bmb.plot.bar(ax=ax2,color='k')
A.plot.bar(ax=ax3,color='k')
ax1.set_ylabel('ismr in m/yr')
ax2.set_ylabel('bmb in Gt/yr')
ax3.set_ylabel('A in 10^3 km^2')
for ax in [ax1,ax2]:
    #ax.set_xticks([])
    ax.set_xticklabels('')
    ax.set_xlabel('')
    ax.xaxis.grid(True)
ax3.xaxis.grid(True)
ax3.set_xlabel('depth bins in m')
plt.savefig(os.path.join(fig_dir,'depth_bins.png'),transparent=True,dpi=300,bbox_inches = "tight")
plt.show()


# ## calculate bmb, ismr and A for shallow and deep for each sector

# In[28]:


sec_depth = {}

for sec_key,sec_mask in log_progress(mask_sector.items(),name='sector'):
    sec_depth[sec_key]={}
    
    for depth_key,depth_mask in mask_depth.items():
        
        mask = sec_mask & depth_mask

        dA = (1/(grd.pm*grd.pn)).where(mask)
        weights = dA/dA.sum()

        ismr2bmb = dA*rhoi*(10**-12)

        sec_depth[sec_key]["A "+depth_key] = dA.sum().values*10**-9
        sec_depth[sec_key]["ismr "+depth_key] = (m.where(mask)*weights).sum().values*s2a
        sec_depth[sec_key]["bmb "+depth_key] = (m.where(mask)*ismr2bmb).sum().values*s2a
        
df_melt = pd.DataFrame.from_dict(sec_depth, orient='index').T


# ## calculate bmb, ismr and A for sizes and sectors

# In[29]:


sec_size = {}

for sec_key,sec_mask in log_progress(mask_sector.items(),name='sector'):
    sec_size[sec_key]={}
    
    for size_key,size_mask in mask_size.items():
        
        mask = sec_mask & size_mask

        dA = (1/(grd.pm*grd.pn)).where(mask)
        weights = dA/dA.sum()

        ismr2bmb = dA*rhoi*(10**-12)

        sec_size[sec_key]["A "+size_key] = dA.sum().values*10**-9
        sec_size[sec_key]["ismr "+size_key] = (m.where(mask)*weights).sum().values*s2a
        sec_size[sec_key]["bmb "+size_key] = (m.where(mask)*ismr2bmb).sum().values*s2a

df_melt = df_melt.append(pd.DataFrame.from_dict(sec_size, orient='index').T)
#df_melt = pd.DataFrame.from_dict(sec_size, orient='index').T


# ## calculate bmb, ismr and A for sizes and sectors for just the front

# In[27]:


mask_front = grd.mask_front==1
sec_size_front = {}

for sec_key,sec_mask in log_progress(mask_sector.items(),name='sector'):
    sec_size_front[sec_key]={}
    
    for size_key,size_mask in mask_size.items():
        
        mask = sec_mask & size_mask & mask_front

        dA = (1/(grd.pm*grd.pn)).where(mask)
        weights = dA/dA.sum()

        ismr2bmb = dA*rhoi*(10**-12)

        sec_size_front[sec_key]["A front "+size_key] = dA.sum().values*10**-9
        sec_size_front[sec_key]["ismr front "+size_key] = (m.where(mask)*weights).sum().values*s2a
        sec_size_front[sec_key]["bmb front "+size_key] = (m.where(mask)*ismr2bmb).sum().values*s2a
        
df_melt = df_melt.append(pd.DataFrame.from_dict(sec_size_front, orient='index').T)


# In[17]:


df_melt = pd.DataFrame.from_dict(sec_size, orient='index').T
df_melt = df_melt.append(pd.DataFrame.from_dict(sec_size_front, orient='index').T)

#df_out_path = os.path.join(int_dir,'melt.csv')
#df_melt.to_csv(df_out_path)


# # plotting

# ## masks

# In[107]:


cma=ocm.deep
colors=(cmap(0.25),cmap(0.5),cmap(0.75))
plt.close()
fig,ax = plt.subplots(figsize=(12,8))
depth_mask.where((grd.mask_rho==1) & (grd.zice <0.0)).plot.contourf(ax=ax,levels=[0,1,2,3],colors=colors)
ax.contourf(grd.mask_rho.where(size==1).values,colors='k',alpha=0.3)
ax.contourf(grd.mask_rho.where(size==2).values,colors='k',alpha=0.5)
plt.show()


# In[112]:


get_ipython().run_line_magic('matplotlib', 'inline')
mask_shallow = (grd.mask_rho==1) & (grd.zice < -200) & (grd.zice >= -400)
mask_deep = (grd.mask_rho==1) & (grd.zice < -400)
levels = [-400,-200]
mask_shelf = (grd.mask_rho == 1) & (grd.zice < 0.0)
plt.close()
fig,ax = plt.subplots(figsize=(8,7))
ax.contourf(grd.mask_rho.where(grd.mask_rho==0).values,colors=(('0.6','0.6','0.6')))
colors=['blue','orange','green','red']
img = size.plot.contourf(ax=ax,levels=[0,1,2,3,4],colors=colors,add_colorbar=False)
ax.contourf(grd.mask_rho.where(mask_shallow).values,colors='k',alpha=0.3)
ax.contourf(grd.mask_rho.where(mask_deep).values,colors='k',alpha=0.5)
#ax.contour(grd.mask_front.where(mask_shelf).values,levels=[0.1],colors='k',linestyles='dashed', alpha=0.75)
#ax.contour(grd.mask_rho.values,levels=['0.01'],colors='k',alpha=1,linewidths=0.5)
#ax.contour(grd.zice.where(grd.mask_rho).values,levels=['-0.01'],colors='k',linestyles='solid',linewidths=0.01)
ax.set_aspect('equal')
ax.axis('off')
#cbaxes = fig.add_axes([0.85, 0.2, 0.02, 0.6])
#cbar = plt.colorbar(img, cax=cbaxes, ticks=[0.5,1.5,2.5,3.5])
#cbar.ax.tick_params(labelsize=18)
#cbar.ax.set_yticklabels(['tiny','small','medium','large'])
plt.tight_layout()
plt.savefig(os.path.join(fig_dir,'size_mask.png'),transparent=True,dpi=300,bbox_to_anchor='tight')
plt.show()


# ## pie charts just bmb

# In[77]:


deep_colors=(ocm.deep(0.25),ocm.deep(0.5),ocm.deep(0.75))

df = df_melt.copy()
df = df.rename(index=str, columns={"Amery/Prydz Bay": "Prydz Bay",
                              "George V Coast/Wilkes subglacial basin":"George V Coast",
                              "Sabrina Coast/Aurora subglacial basin":"Sabrina Coast"})

plt.close()
fig,ax = plt.subplots(figsize=(10,10))
radius = df['Total Antarctica'].loc['bmb all'] / df['Bellingshausen Sea'].loc['bmb all']/4
img = df['Total Antarctica'].loc[['bmb shallow','bmb mid','bmb deep']].plot.pie(ax=ax,labels=None,radius=radius,colors=deep_colors,
                                                                                      wedgeprops={"edgecolor":"k",'linewidth': 5, 'antialiased': True})
ax.set_title('Total Antarctica',fontsize=40)
ax.set_aspect('equal')
ax.set_ylabel('')
plt.tight_layout()
out_path= os.path.join(fig_dir,'tmp','depth_pie_Total_Antarctica.png')

plt.savefig(out_path,format='png',transparent=True,bbox_inches = "tight")

for sec_name, sec_data in df.drop(columns='Total Antarctica').items():
    fig,ax = plt.subplots(figsize=(10,10))
    radius = sec_data.loc['bmb all'] / df['Bellingshausen Sea'].loc['bmb all']
    sec_data.loc[['bmb shallow','bmb mid','bmb deep']].plot.pie(ax=ax,radius=radius,labels=None,colors=deep_colors,
                                                                            wedgeprops={"edgecolor":"k",'linewidth': 5, 'antialiased': True})
    ax.set_aspect('equal')
    ax.set_title(sec_name,fontsize=40)
    ax.set_ylabel('')
    out_path= os.path.join(fig_dir,'tmp','depth_pie_'+sec_name.replace('/','_')+'.png')
    plt.tight_layout()
    plt.savefig(out_path,format='png',transparent=True,bbox_inches = "tight")
    if sec_name == 'Prydz Bay':
        fig,ax = plt.subplots(figsize=(10,10))
        radius = sec_data.loc['bmb all'] / df['Bellingshausen Sea'].loc['bmb all'] * 5.0
        sec_data.loc[['bmb shallow','bmb mid','bmb deep']].plot.pie(ax=ax,radius=radius,labels=None,colors=deep_colors,
                                                                            wedgeprops={"edgecolor":"k",'linewidth': 5, 'antialiased': True})
        ax.set_aspect('equal')
        ax.set_title(sec_name,fontsize=40)
        ax.set_ylabel('')
        out_path= os.path.join(fig_dir,'tmp','depth_pie_'+sec_name.replace('/','_')+'_x5.png')
        plt.tight_layout()
        plt.savefig(out_path,format='png',transparent=True,bbox_inches = "tight")


# ## Table mass loss, area and ismr for diffrent depths

# In[13]:


shelves_melt = shelves.copy()


# In[22]:


shelves_melt = {}
for shelf_key,shelf_data in log_progress(shelves.items(),name='ice shelf'):
    shelf_mask=shelf_data['mask']
    shelves_melt[shelf_key]={}
    
    for depth_key,depth_mask in mask_depth.items():
        
        mask = shelf_mask & depth_mask

        dA = (1/(grd.pm*grd.pn)).where(mask)
        weights = dA/dA.sum()

        ismr2bmb = dA*rhoi*(10**-12)

        shelves_melt[shelf_key]["A "+depth_key] = dA.sum().values*10**-9
        shelves_melt[shelf_key]["ismr "+depth_key] = (m.where(mask)*weights).sum().values*s2a
        shelves_melt[shelf_key]["bmb "+depth_key] = (m.where(mask)*ismr2bmb).sum().values*s2a
        
shelves_df = pd.DataFrame.from_dict(shelves, orient='index').T


# In[23]:


df_shelves =  pd.DataFrame.from_dict(shelves_melt,orient='index')


# In[24]:


df_shelves


# In[30]:


df_sectors = df_melt.copy()


# In[31]:


def wavg(group, avg_name, weight_name):
    d = group[avg_name]
    w = group[weight_name]
    try:
        return (d * w).sum() / w.sum()
    except ZeroDivisionError:
        return d.mean()
    
def combine_shelves(df,name):
    
    new_data = {}
    
    for depth in ['all','<200','<400','>400']:
        key = 'ismr '+depth
        new_data[key]=wavg(df,key,'A '+depth)
        for quant in ['bmb','A']:
            key = quant+' '+depth
            new_data[key]=df[key].sum()
        
    new_dict = {name:{'bmb all':new_data['bmb all'],'ismr all':new_data['ismr all'],'A all':new_data['A all'],
                      'bmb <200':new_data['bmb <200'],'ismr <200':new_data['ismr <200'],'A <200':new_data['A <200'],
                      'bmb <400':new_data['bmb <400'],'ismr <400':new_data['ismr <400'],'A <400':new_data['A <400'],
                      'bmb >400':new_data['bmb >400'],'ismr >400':new_data['ismr >400'],'A >400':new_data['A >400']}}
                      
    new_df = pd.DataFrame.from_dict(new_dict, orient='index')
    
    return new_df


# In[32]:


ant_bounds_comb = combine_shelves(df_shelves,'Ant_bounds all')
ant_bounds_comb


# In[45]:


df_melt[0:9].T


# In[52]:


out = pd.concat([df_shelves,ant_bounds_comb,df_melt[0:9].T.sort_index(by='bmb all')])


# In[53]:


for quant in ['bmb','A']:
    for depth in ['<200','<400','>400']:    
        key = quant+' '+depth 
        out[key] = out[key]/out[quant+' all']*100
        out[quant] = out[quant+' all']
    out=out.drop(quant+' all',axis=1)
    
out['ismr'] = out['ismr all']


# In[ ]:





# In[54]:


out = out[['ismr','A','bmb','ismr <200','A <200','bmb <200','ismr <400','A <400','bmb <400']].dropna()
pd.options.display.float_format = '{:,.2f}'.format
out.to_csv(os.path.join(proj_dir,'reports','tables','2km_mass_loss_bydepth.csv'))
out


# In[27]:


out['bmb shallow']


# ## histogram plots for A, bmb and ismr

# In[25]:


A_max = df.drop(columns='Total Antarctica').T['A all'].max()
bmb_max = df.drop(columns='Total Antarctica').T['bmb all'].max()
ismr_max = df.drop(columns='Total Antarctica').T['ismr all'].max()


# In[41]:


plt.close()
for sec_name,sec_data in df.drop(columns='Total Antarctica').items():
    fig,axes = plt.subplots(ncols=3)
    fig.suptitle(sec_name)
    for ax,key in zip(axes.flatten(),['A','bmb','ismr']):
        sec_data.loc[[key+' all']].plot.bar(ax=ax,use_index=False,color='k')
        ax.set_ylim(0,df.drop(columns='Total Antarctica').T[key+' all'].max())
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(key)
    #plt.tight_layout()    
    plt.savefig(os.path.join(fig_dir,'tmp','all_'+sec_key.replace('/','_'))+'.png',format='png',transparent=True)
    plt.show()


# In[332]:


plt.close()
fig = plt.figure(figsize=(15,10))
df.loc[['bmb tiny','bmb small','bmb medium','bmb large','bmb all']].T.plot.bar(legend=True)
plt.legend(['Tiny (Not in Antbound)','Small (0-5 *1000 km^2)','Medium (5-45 *1000 km^2)','Large (>45 *1000 km^2)'],fontsize=35,
          bbox_to_anchor=(1,1))
plt.savefig(os.path.join(fig_dir,'tmp','pie_legend.png'),transparent=True,bbox_inches='tight')
plt.show()


# ## continental shelf temps vs ismr

# ## calculate average on shelf temps (shallow, deep, cavity) for sectors

# In[7]:


#%%writefile ../../src/tools/make_grd_dV.py
# equip roms grid with depths and volumes
from tools.cartesian_grid_3d import cartesian_grid_3d
import numpy as np
import xarray as xr

def make_grd_dV(grd,zeta):

    lon_u = grd.lon_u.values
    lat_u = grd.lat_u.values
    lon_v = grd.lon_v.values
    lat_v = grd.lat_v.values
    h = grd.h.values
    zice = grd.zice.values
    theta_s = 7#zeta.theta_s.values
    theta_b = 8#zeta.theta_b.values
    hc = 250#temp.hc.values
    N = 31

    dx,dy,dz,z = cartesian_grid_3d(lon_u, lat_u, lon_v, lat_v, h, zice, theta_s, theta_b, hc, N, zeta.values)

    grd['dx'] = xr.DataArray(dx,dims=['s_rho','eta_rho','xi_rho'])
    grd['dx'] = grd.dx.where(grd.mask_rho == 1)

    grd['dy'] = xr.DataArray(dy,dims=['s_rho','eta_rho','xi_rho'])
    grd['dy'] = grd.dy.where(grd.mask_rho == 1)

    grd['dz'] = xr.DataArray(dz,dims=['s_rho','eta_rho','xi_rho'])
    grd['dz'] = grd.dz.where(grd.mask_rho == 1)

    grd['z'] = xr.DataArray(z,dims=['s_rho','eta_rho','xi_rho'])
    grd['z'] = grd.z.where(grd.mask_rho == 1)

    dV = grd.dx * grd.dy * grd.dz
    grd['dV'] = dV

    return grd


# In[6]:


grd=make_grd_dV(grd)
grd_out_path = os.path.join(int_dir,'grd4_dV.nc')
grd.to_netcdf(grd_out_path)


# In[77]:


mask_cavity = (grd.mask_rho == 1) & (grd.zice < 0.0)
mask =mask_shelf_sector['Amery/Prydz Bay']

plt.close()
mask.plot()
plt.show()


# In[115]:


sec_depth_temp = {}

mask_cavity =  (grd.mask_rho == 1) & (grd.zice < 0.0)
mask_shelf = mask_shelf_sector['Total Antarctica']

for sec_name,sec_mask in log_progress(mask_sector.items()):
    sec_depth_temp[sec_name]={}
    
    for depth_name,depth_mask in zip(['shallow'],[(mask_shelf | mask_cavity) & (grd.z >= -200)]):
    #zip(['all','cavity','shelf','100','200','350','500'],
                                 #    [mask_cavity | mask_shelf,mask_cavity,mask_shelf,
                                 #     (mask_shelf == 1) & (grd.z < -100),
                                 #     (mask_shelf == 1) & (grd.z < -200),
                                 #    (mask_shelf == 1) & (grd.z < -350),
                                 #    (mask_shelf == 1) & (grd.z < -500)]): 
    
        mask = sec_mask & depth_mask
    
        dV = grd.dV.where(mask)
        weights = dV/dV.sum()
    
        sec_depth_temp[sec_name][depth_name] = (temp.where(mask)*weights).sum().values*1.0


# In[119]:


for sec_name,sec_mask in log_progress(mask_shelf_sector.items()):
    
    mask = (sec_mask==1) & (grd.zice < 0.0) & (grd.mask_rho==1)
    
    sec_depth_temp[sec_name]['avg draft'] = grd.zice.where(mask).mean().values*1.0


# In[120]:


df_temp = pd.DataFrame.from_dict(sec_depth_temp, orient='index').T

df_out_path = os.path.join(int_dir,'temp.csv')
df_temp.to_csv(df_out_path)


# In[121]:


df = pd.concat([df_melt,df_temp])


# In[124]:


df_plt = df.T.drop('Total Antarctica').sort_values('ismr all',ascending=False)
plt.close()
fig,ax = plt.subplots(figsize=(10,4))
df_plt.plot(y='ismr front all',ax=ax)
df_plt.plot(y='shallow',linestyle='dashed',ax=ax,xticks=range(8),rot=90,secondary_y=True)
#df_plt.plot(y='avg draft',ax=ax2,xticks=range(8),rot=90)
#df_plt.plot(y='ismr front all',ax=ax,xticks=range(8),rot=90,secondary_y=True)

#plt.savefig(os.path.join(fig_dir,'tmp','ismr_vs_temps.png'),bbox_inches='tight')
plt.show()


# In[39]:


df_plt = df.T.drop('Total Antarctica').sort_values('ismr all',ascending=False)
plt.close()
fig,axes = plt.subplots(ncols=2,figsize=(10,4))
ax1,ax2 = axes.flatten()
df_plt.plot(y=['ismr all','ismr front all'],ax=ax1)
df_plt.plot(y=['shallow','deep','cavity'],linestyle='dashed',ax=ax1,xticks=range(8),rot=90,secondary_y=True)
df_plt.plot(y='avg draft',ax=ax2,xticks=range(8),rot=90)
#df_plt.plot(y='ismr front all',ax=ax,xticks=range(8),rot=90,secondary_y=True)

plt.savefig(os.path.join(fig_dir,'tmp','ismr_vs_temps.png'),bbox_inches='tight')
plt.show()


# In[52]:


plt.close()
fig,ax = plt.subplots()
df.T.plot.bar(y=['avg draft'],colors='k',ax=ax)
#df.T.plot.bar(y=['avg draft'],secondary_y=True,colors='k',ax=ax)

plt.show()


# In[376]:


def make_sector(mask_sector,name_sector):
    
    sector = {}
    
    mask = mask_sector & (grd.mask_rho==1) & (grd.zice<0)

    dA = (1/(grd.pm*grd.pn)).where(mask)
    weights = dA/dA.sum()

    ismr2bmb = dA*rhoi*(10**-12)

    sector["A total"] = dA.sum().values*10**-9
    sector["ismr total"] = (m.where(mask)*weights).sum().values*s2a
    sector["bmb total"] = (m.where(mask)*ismr2bmb).sum().values*s2a
    
    for k,v in mask_size.items():
        
        mask = mask_sector & (grd.mask_rho==1) & (grd.zice<0) & (v==1)

        dA = (1/(grd.pm*grd.pn)).where(mask)
        weights = dA/dA.sum()

        ismr2bmb = dA*rhoi*(10**-12)

        sector["A "+k] = dA.sum().values*10**-9
        sector["ismr "+k] = (m.where(mask)*weights).sum().values*s2a
        sector["bmb "+k] = (m.where(mask)*ismr2bmb).sum().values*s2a
    
    sector_df = pd.DataFrame.from_dict({name_sector:sector}, orient='index').dropna()
    
    return sector_df


# In[389]:


for k,v in mask_sectors.items():  
    all_ice_df = all_ice_df.append(make_sector(v,k))


# In[390]:


all_ice_df


# In[324]:


mask_Ronne = make_mask_lonlat(grd,-85,-49.5,-84.50,-74.7)
mask_Filchner = make_mask_lonlat(grd,-49.5,-27.66,-83.5,-77.73)
mask_FRIS = mask_Filchner | mask_Ronne

mask_Amery = make_mask_lonlat(grd,65,74.3,-73.7,-68.3) 

mask_RossIS_West = make_mask_lonlat(grd,-180,-146.70,-86,-77.8)
mask_RossIS_East = make_mask_lonlat(grd,158.3,180,-84.5,-77)
mask_RossIS = mask_RossIS_West | mask_RossIS_East

mask_Larsen = make_mask_lonlat(grd,-65.5,-60,-69.3,-66.1)

vostock = (grd.lat_rho<-75) & (grd.lat_rho>-80) & (grd.lon_rho>95) & (grd.lon_rho<115)

mask_BA = make_mask_lonlat(grd,-141,-65.5,-90,-60) & ~mask_FRIS
mask_Ross = make_mask_lonlat(grd,160,-141,-90,-60)
mask_EI = make_mask_lonlat(grd,94.2,160,-90,-60) & ~vostock
mask_WI = make_mask_lonlat(grd,7.6,94.2,-90,-60)
mask_Weddell = make_mask_lonlat(grd,-65.5,7.6,-90,-60) | mask_FRIS


# In[230]:


all_ice = {}

ice_shelf = (grd.mask_rho==1) & (grd.zice<0.0) 
vostock = (grd.lat_rho<-75) & (grd.lat_rho>-80) & (grd.lon_rho>95) & (grd.lon_rho<115)

mask = ice_shelf & ~vostock

dA = (1/(grd.pm*grd.pn)).where(mask)
weights = dA/dA.sum()

ismr2bmb = dA*rhoi*(10**-12)

all_ice["A"] = dA.sum().values*10**-9
all_ice["ismr"] = (m.where(mask)*weights).sum().values*s2a
all_ice["bmb"] = (m.where(mask)*ismr2bmb).sum().values*s2a

mask = mask & (grd.mask_front == 0)

dA = (1/(grd.pm*grd.pn)).where(mask)
weights = dA/dA.sum()

ismr2bmb = dA*rhoi*(10**-12)

all_ice["A_nf"] = dA.sum().values*10**-9
all_ice["ismr_nf"] = (m.where(mask)*weights).sum().values*s2a
all_ice["bmb_nf"] = (m.where(mask)*ismr2bmb).sum().values*s2a

mask = (grd.h < 1000.0) & (grd.zice == 0.0) & (grd.mask_rho == 1)
    
dA = (1/(grd.pm*grd.pn)).where(mask)
weights = dA/dA.sum()
    
all_ice["theta"] = (temp.where(mask)*weights).sum().values

all_ice = pd.DataFrame.from_dict({'All ice (Bedmap2)':all_ice}, orient='index',
                            columns=['bmb','bmb_nf','ismr','ismr_nf','A','A_nf','theta']).dropna()
all_ice


# # Ice by ocean sector

# In[226]:


def make_mask_lonlat(grd,lon_min,lon_max,lat_min,lat_max):

    if lon_min<lon_max:
        mask_lonlat = (grd.lon_rho>lon_min) & (grd.lon_rho<=lon_max) & (grd.lat_rho>lat_min) & (grd.lat_rho<=lat_max)
    else:
        mask_lonlat = ((grd.lon_rho>lon_min) | (grd.lon_rho<=lon_max)) & (grd.lat_rho>lat_min) & (grd.lat_rho<=lat_max)

    
    return mask_lonlat  


# In[278]:


mask_Ronne = make_mask_lonlat(grd,-85,-49.5,-84.50,-74.7)
mask_Filchner = make_mask_lonlat(grd,-49.5,-27.66,-83.5,-77.73)
mask_FRIS = mask_Filchner | mask_Ronne

mask_Amery = make_mask_lonlat(grd,65,74.3,-73.7,-68.3) 

mask_RossIS_West = make_mask_lonlat(grd,-180,-146.70,-86,-77.8)
mask_RossIS_East = make_mask_lonlat(grd,158.3,180,-84.5,-77)
mask_RossIS = mask_RossIS_West | mask_RossIS_East

mask_Larsen = make_mask_lonlat(grd,-65.5,-60,-69.3,-66.1)

vostock = (grd.lat_rho<-75) & (grd.lat_rho>-80) & (grd.lon_rho>95) & (grd.lon_rho<115)

mask_BA = make_mask_lonlat(grd,-141,-65.5,-90,-60) & ~mask_FRIS
mask_Ross = make_mask_lonlat(grd,160,-141,-90,-60)
mask_EI = make_mask_lonlat(grd,94.2,160,-90,-60) & ~vostock
mask_WI = make_mask_lonlat(grd,7.6,94.2,-90,-60)
mask_Weddell = make_mask_lonlat(grd,-65.5,7.6,-90,-60) | mask_FRIS


# In[279]:


get_ipython().run_line_magic('matplotlib', 'inline')



mask_sec = mask_Weddell

mask_ice =  (grd.mask_rho==1) & (grd.zice<0)
mask_temp =  (grd.h < 1000.0) & (grd.mask_rho == 1) & (grd.zice == 0.0)

plt.close()
fig,ax = plt.subplots(figsize=(15,10))
temp.where(mask_temp).plot(ax=ax,vmin=-2,vmax=1)
(m.where(mask_ice)*365*24*3600).plot(ax=ax,vmin=-2,vmax=2)

for mask_sec in [mask_BA,mask_Ross,mask_EI,mask_WI,mask_Weddell]: 
    ax.contour(mask_sec)

ax.set_aspect('equal')
plt.show()


# In[244]:


def make_row(mask_ice_sec,mask_temp_sec,name):
    
    row = {}
    
    mask = mask_ice_sec & (grd.mask_rho==1) & (grd.zice<0)

    dA = (1/(grd.pm*grd.pn)).where(mask)
    weights = dA/dA.sum()

    ismr2bmb = dA*rhoi*(10**-12)

    row["A"] = dA.sum().values*10**-9
    row["ismr"] = (m.where(mask)*weights).sum().values*s2a
    row["bmb"] = (m.where(mask)*ismr2bmb).sum().values*s2a
    
    mask = mask & (grd.mask_front==0)
    
    dA = (1/(grd.pm*grd.pn)).where(mask)
    weights = dA/dA.sum()

    ismr2bmb = dA*rhoi*(10**-12)
    
    row["A_nf"] = dA.sum().values*10**-9
    row["ismr_nf"] = (m.where(mask)*weights).sum().values*s2a
    row["bmb_nf"] = (m.where(mask)*ismr2bmb).sum().values*s2a
    
    mask = mask_temp_sec & (grd.h < 1000.0) & (grd.zice == 0.0) & (grd.mask_rho == 1)
    
    dA = (1/(grd.pm*grd.pn)).where(mask)
    weights = dA/dA.sum()
    
    row["theta"] = (temp.where(mask)*weights).sum().values
    
    row = pd.DataFrame.from_dict({name:row}, orient='index',
                                columns=['bmb','bmb_nf','ismr','ismr_nf','A','A_nf','theta']).dropna()
    return row


# In[247]:


regime = make_row(mask_BA,mask_BA,'BA all')
regime = regime.append(make_row(mask_lc,(mask_Ross | mask_WI | mask_Weddell),'large cold'))
regime = regime.append(make_row(mask_sc,~mask_BA,'small cold'))
regime


# In[271]:


sectors = make_row(mask_BA,mask_BA,'BA Seas')

for mask_ice,mask_temp,name in log_progress(zip([mask_RossIS,(mask_Ross & ~mask_RossIS),
                                                 mask_EI,
                                                 mask_Amery, (mask_WI & ~mask_Amery),
                                                 mask_FRIS,mask_Larsen,(mask_Weddell & ~mask_FRIS)],
                                                [mask_Ross,mask_Ross,
                                                 mask_EI,
                                                 mask_WI,mask_WI,
                                                 mask_Weddell,mask_Weddell,mask_Weddell],
                                                ['Ross','Ross sea small',
                                                 'East Indian sea',
                                                 'Amery','West Indian sea small',
                                                 'FRIS','LarsenC','Weddell Sea small']),every=3):
    sectors = sectors.append(make_row(mask_ice,mask_temp,name))
sectors


# In[272]:


pd.options.display.float_format = '{:,.2f}'.format
out_dir = os.path.join(pro_dir,'2km_mass_loss_ocean_sectors.csv')
sectors.to_csv(out_dir)


# In[268]:


out = pd.concat([all_ice,regime])
out['bmb/bmb_total'] = out['bmb']/1209
out['bmb_nf/bmb_nf_total'] = out['bmb_nf']/829
out = out[['bmb','bmb/bmb_total','ismr','A','bmb_nf','bmb_nf/bmb_nf_total','ismr_nf','A_nf']]
out


# In[269]:


pd.options.display.float_format = '{:,.2f}'.format
out_dir = os.path.join(pro_dir,'2km_mass_loss_size_regime.csv')
out.to_csv(out_dir)


# ## All ice shelves from antbound

# In[168]:


#optional plot ice shelf masks and names
get_ipython().run_line_magic('matplotlib', 'notebook')

plt.close()
fig,ax = plt.subplots(figsize=(10,7))
grd.zice.where((grd.mask_rho==1)&(grd.zice<0)).plot(ax=ax,alpha=0.2,add_colorbar=False)

for name,data in log_progress(shelves.items()):
    
    grd.zice.where(data['mask']).plot(ax=ax,add_colorbar=False,alpha=0.5)
    ax.text(data['xi_mid'],data['eta_mid'],name,alpha=0.8)
    
ax.set_aspect('equal')
plt.tight_layout()

plt.show()


# In[334]:


df_out_path = os.path.join(int_dir,'melt.csv')
df_melt.to_csv(df_out_path)#avg = xr.open_dataset(avg_path)

#m = avg.m.mean('ocean_time')

#s2a = 3600*24*365.25
#rhoi = 916

mask_ice = (grd.mask_rho==1) & (grd.zice<0)

for name,data in log_progress(shelves.items(),name='Ice shelf'): 

    mask = data['mask']

    dA = (1/(grd.pm*grd.pn)).where(mask)
    weights = dA/dA.sum()

    dA_l = dA.where(m > 0.0)
    weights_l = dA_l/dA.sum()

    dA_g = dA.where(m < 0.0)
    weights_g = dA_g/dA.sum()

    ismr2bmb = dA*rhoi*(10**-12)

    shelves[name]["A"] = dA.sum().values*10**-9
    shelves[name]["ismr"] = (m.where(mask)*weights).sum().values*s2a
    shelves[name]["ismr_l"] = (m.where(mask & (m > 0.0))*weights_l).sum().values*s2a
    shelves[name]["ismr_g"] = (m.where(mask & (m < 0.0))*weights_g).sum().values*s2a
    shelves[name]["bmb"] = (m.where(mask)*ismr2bmb).sum().values*s2a
    shelves[name]["bml"] = (m.where(mask & (m > 0.0))*ismr2bmb).sum().values*s2a
    shelves[name]["bmg"] = (m.where(mask & (m < 0.0))*ismr2bmb).sum().values*s2a
    shelves[name]['ismr_max'] = m.where(mask).max().values*s2a
    shelves[name]['ismr_min'] = m.where(mask).min().values*s2a


# In[336]:


df = pd.DataFrame.from_dict(shelves, orient='index',
                            columns=['bmb','ismr','A','bml','bmg','ismr_l','ismr_g','ismr_max','ismr_min']).dropna()


# In[337]:


pd.options.display.float_format = '{:,.2f}'.format
df.to_csv(os.path.join(int_dir,'2km _mass_loss_ant_bounds_all.csv'))
df


# In[29]:


def wavg(group, avg_name, weight_name):
    d = group[avg_name]
    w = group[weight_name]
    try:
        return (d * w).sum() / w.sum()
    except ZeroDivisionError:
        return d.mean()
    
def combine_shelves(df,name):
    
    new_data = {}
    
    for key in ['bmb','bml','bmg','A']:
        new_data[key]=df[key].sum()
    for key in ['ismr','ismr_l','ismr_g']:
        new_data[key]=wavg(df,key,'A')
        
    new_data['ismr_max']=df['ismr_max'].max()
    new_data['ismr_min']=df['ismr_min'].min()
    
    new_dict = {name:{'bmb':new_data['bmb'],'ismr':new_data['ismr'],'A':new_data['A'],
                      'bml':new_data['bml'],'bmg':new_data['bmg'],
                      'ismr_l':new_data['ismr_l'],'ismr_g':new_data['ismr_g'],
                      'ismr_max':new_data['ismr_max'],'ismr_min':new_data['ismr_min']}}
    new_df = pd.DataFrame.from_dict(new_dict, orient='index',
                                    columns=['bmb','ismr','A','bml','bmg',
                                             'ismr_l','ismr_g','ismr_max','ismr_min'])
    
    return new_df


# In[355]:


ant_bounds_comb = combine_shelves(df,'Ant_bounds all')
ant_bounds_comb


# # following Rignot 2013

# In[340]:


def select_names(df,names):
    df_select = df[df.index.str.contains('bla')]
    for string in names:
        df_select = df_select.append(df[df.index.str.contains(string)])
    return df_select


# In[341]:


AP_names = ['LarsenG','LarsenF','LarsenE','LarsenD','LarsenC','LarsenB','Wordie','Wilkins','Bach','George VI','Stange']
select_names(df,AP_names)


# In[342]:


AP = pd.concat([select_names(df,['LarsenG','LarsenF','LarsenE']),
              combine_shelves(select_names(df,['LarsenD']),'Larsen D'),
              select_names(df,['LarsenC']),
              combine_shelves(select_names(df,['Wordie']),'Wordie'),
              select_names(df,['Wilkins','Bach','George VI','Stange'])])
AP


# In[343]:


AP_comb = combine_shelves(AP,'AP')
AP_comb


# In[344]:


WAIS_names = ['Ronne','Ferrigno','Venable','Abbot','Cosgrove','Pine Island','Thwaites','Crosson','Dotson','Getz',
              'Land','Nickerson','Sulzberger','Swinburne','Withrow','Ross West']
select_names(df,WAIS_names)


# In[345]:


WAIS = pd.concat([select_names(df,['Ronne','Ferrigno','Venable']),
              combine_shelves(select_names(df,['Abbot']),'Abbot'),
              select_names(df,['Cosgrove','Pine Island','Thwaites','Crosson','Dotson']),
              combine_shelves(select_names(df,['Getz']),'Getz'),
              select_names(df,['Nickerson','Sulzberger','Swinburne','Withrow','Ross West'])])
WAIS


# In[346]:


WAIS_comb = combine_shelves(WAIS,'WAIS')
WAIS_comb


# In[347]:


EAIS_names = ['Ross East','Drygalski','Nansen','Aviator','Mariner','Lillie','Rennick','Cook','Ninnis','Mertz',
              'Dibble','Holmes','Moscow University','Totten','Vincennes','Conger','Tracy','Shackleton','West',
              'Publications','Amery','Wilma','Edward','Rayner','Shirase','Prince Harald','Baudouin','Borchgrevink',
              'Lazarev','Nivl','Vigrid','Fimbul','Jelbart','Atka','Ekstrom','Quar','Riiser','Brunt','Filchner']
select_names(df,EAIS_names)


# In[348]:


EAIS = select_names(df,EAIS_names).drop('Ross West')
EAIS


# In[349]:


EAIS_comb = combine_shelves(EAIS,'EAIS')
EAIS_comb


# In[352]:


Rignot_comb = combine_shelves(pd.concat([AP_comb,WAIS_comb,EAIS_comb]),'Rignot all')
Rignot_comb


# In[356]:


df_out = pd.concat([AP,AP_comb,
           WAIS,WAIS_comb,
           EAIS,EAIS_comb,
           Rignot_comb,
           ant_bounds_comb,
           all_ice])
df_out


# In[357]:


csv_out_path = os.path.join(pro_dir,'2km_mass_loss.csv')
df_out.to_csv(csv_out_path)


# # Major ice shelves for comparison

# In[245]:


names = ['George','Abbot','Pine','Getz','Fimbul','Jelbart','Brunt','Riiser','Filchner','Ronne','LarsenC','Ross',
                   'Totten','Moscow','Shackleton','Amery']
major = select_names(df,names)     


# In[246]:


def merge_shelves(df,name1,name2,name):
    
    p1 = df.loc[name1]
    p2 = df.loc[name2]
    
    new_data = {}
    for key in ['bmb','bml','bmg','A']:
        new_data[key]=p1[key]+p2[key]
    for key in ['ismr','ismr_l','ismr_g']:
        new_data[key]=(p1[key]*p1['A']+p2[key]*p2['A'])/(p1['A']+p2['A'])
        
    new_data['ismr_max']=np.max([p1['ismr_max'],p2['ismr_max']])
    new_data['ismr_min']=np.max([p1['ismr_min'],p2['ismr_min']])
    
    new_dict = {name:{'bmb':new_data['bmb'],'ismr':new_data['ismr'],'A':new_data['A'],
                      'bml':new_data['bml'],'bmg':new_data['bmg'],
                      'ismr_l':new_data['ismr_l'],'ismr_g':new_data['ismr_g'],
                      'ismr_max':new_data['ismr_max'],'ismr_min':new_data['ismr_min']}}
    new_df = pd.DataFrame.from_dict(new_dict, orient='index',columns=['bmb','ismr','A','bml','bmg','ismr_l','ismr_g','ismr_max','ismr_min'])
    df = df.drop([name1,name2])
    df = df.append(new_df)
    
    return df


# In[247]:


major = merge_shelves(major,'Abbot','Abbot 1','Abbot')
for key in ['2','3','4','5','6']:
    major = merge_shelves(major,'Abbot','Abbot '+key,'Abbot')
major


# In[248]:


major = merge_shelves(major,'Filchner','Ronne','Filchner-Ronne')
major = merge_shelves(major,'Getz','Getz 1','Getz')
major = merge_shelves(major,'Fimbul','Jelbart','Fimbul + Jelbart')
major = merge_shelves(major,'Brunt Stancomb','Riiser-Larsen','Brunt + Riiser-Larsen')
major = merge_shelves(major,'Ross East','Ross West','Ross')
major = merge_shelves(major,'Totten','Moscow University','Totten + Moscow Uni')


# In[249]:


major


# In[358]:


csv_select_pro_path = os.path.join(pro_dir,'2km_mass_loss_major.csv')
major.to_csv(csv_select_pro_path)


# In[201]:


major_comb = combine_shelves(major,'major')
major_comb


# # Major Ice shelves after ocean and size

# In[203]:


names = ['George','Wilkins','Bach','Stange','Ferrigno',
        'Venable','Abbot','Cosgrove','Pine','Thwaits','Crosson',
         'Dotson','Getz','Land']
warm_AB = select_names(df,names)
warm_AB


# In[230]:


warm_AB_comb = combine_shelves(warm_AB,'Amunden-Bellinghausen')
warm_AB_comb


# In[205]:


names = ['Filchner','Ronne','Ross','LarsenC','Amery']
large_cold = select_names(df,names)
large_cold


# In[207]:


large_cold_comb = combine_shelves(large_cold,'large cold')
large_cold_comb


# In[180]:


names = ['Nickerson','Sulzberg','Swinburne','Withrow','Drygal','Nansen',
         'Aviator','Mariner']
small_RS = select_names(df,names)
small_RS


# In[210]:


small_RS_comb = combine_shelves(small_RS,'Small Ross Sea')
small_RS_comb


# In[185]:


names = ['Lillie','Rennick','Cook','Mertz','Dibble',
         'Holmes','Moscow','Totten','Vincennes','Conger','Tremenchus',
        'Shackleton']
small_EI = select_names(df,names)
small_EI


# In[211]:


small_EI_comb = combine_shelves(small_EI,'small East Indien')
small_EI_comb


# In[190]:


names = ['Vigrid','Nivl','Lazarev','Borch','Bau','Prince','Shirase',
         'Rayner','Edward','Wilma','Publications','West']
small_WI = select_names(df,names).drop('Ross West')
small_WI


# In[212]:


small_WI_comb = combine_shelves(small_WI,'Small West Indian')
small_WI_comb


# In[193]:


names = ['LarsenB','LarsenD','LarsenE','LarsenF','LarsenG',
         'Stancomb','Riiser','Quar','Ekstr','Atka','Jelbart','Fimbul']
small_WS = select_names(df,names)
small_WS


# In[213]:


small_WS_comb = combine_shelves(small_WS,'Small Weddell Sea')
small_WS_comb


# In[236]:


comb_all = combine_shelves(pd.concat([warm_AB_comb,large_cold_comb,small_EI_comb,small_WI_comb,small_WS_comb]),'combined all')
comb_small = combine_shelves(pd.concat([small_EI_comb,small_WI_comb,small_WS_comb]),'combined small')


# In[359]:


regime_df = pd.concat([warm_AB_comb,large_cold_comb,small_EI_comb,small_WI_comb,small_WS_comb,comb_small,comb_all])
out_path = os.path.join(pro_dir,'2km_mass_loss_regime.csv')
regime_df.to_csv(out_path)


# In[104]:


out_path = os.path.join(fig_dir,'bmb2_select.png')

plt.close()
fig,ax = plt.subplots(figsize=(15,30))
df_select.plot(y=['bmb','bml','bmg'],kind='barh',ax=ax,fontsize=24,color=['black','gray','gray'])
ax.legend(markerscale=1,fontsize=24,loc=4)
ax.set_title('Ice shelf basal mass balance',fontsize=24)
ax.set_xlabel('Mass change rate in Gt/yr',fontsize=24)
ax.grid(axis='x')
#plt.savefig(out_path,dpi=300,format='png', bbox_inches='tight')
plt.show()


# In[31]:


out_path = os.path.join(fig_dir,'ismr_all_2.png')

plt.close()
fig,ax = plt.subplots(figsize=(15,30))
df.plot(y=['ismr','ismr_l','ismr_g'],kind='barh',ax=ax,fontsize=24,color=['black','gray','gray'])
ax.legend(markerscale=1,fontsize=24,loc=4)
ax.set_title('Ice shelf basal mass balance',fontsize=24)
ax.set_xlabel('Area average melt rate in m/yr',fontsize=24)
ax.grid(axis='x')
plt.savefig(out_path,dpi=300,format='png', bbox_inches='tight')
plt.show()


# In[ ]:





# In[ ]:




def make_mass_loss_table(m,grd,
                         google_sheet_url='https://docs.google.com/spreadsheets/d/1BpI6iRB569kF7TxdiFT3ColNVH9zwHgNwq49Dat7Ok4/edit?usp=sharing',
                         ,include_no_front=True):
    
    csv_export_url = google_sheet_url.replace('edit?', 'export?gid=0&format=csv&')
    IS = pd.read_csv(csv_export_url)

    s2a = 3600*24*365
    rhoi = 916

    mask_ice = (grd.mask_rho==1) & (grd.zice<0)

    for idx,row in log_progress(IS.iloc[idx_start:idx_end].iterrows(),name='Ice shelf',every=1): 

        if row.lon_min<row.lon_max:
            mask_coord = (grd.lon_rho>row.lon_min) & (grd.lon_rho<=row.lon_max) & (grd.lat_rho>row.lat_min) & (grd.lat_rho<=row.lat_max)
        else:
            mask_coord = ((grd.lon_rho>row.lon_min) | (grd.lon_rho<=row.lon_max)) & (grd.lat_rho>row.lat_min) & (grd.lat_rho<=row.lat_max)

        mask_shelf = mask_ice & mask_coord
        
        dA = (1/(grd.pm*grd.pn)).where(mask_shelf)
        weights = dA/dA.sum()
        
        dA_l = dA.where(m > 0.0)
        weights_l = dA_l/dA_l.sum()
        
        dA_g = dA.where(m < 0.0)
        weights_g = dA_g/dA_g.sum()
        
        ismr2bmb = dA*rhoi*(10**-12)

        IS.at[idx,"A"] = float(dA.sum().values)*10**-9
        IS.at[idx,"ismr"] = float((m.where(mask_shelf)*weights*s2a).sum().values)
        IS.at[idx,"ismr_l"] = float((m.where(mask_shelf & (m > 0.0))*weights_l*s2a).sum().values)
        IS.at[idx,"ismr_g"] = float((m.where(mask_shelf & (m < 0.0))*weights_g*s2a).sum().values)
        IS.at[idx,"bmb"] = float((m.where(mask_shelf)*ismr2bmb*s2a).sum().values)
        IS.at[idx,"bml"] = float((m.where(mask_shelf & (m > 0.0))*ismr2bmb*s2a).sum().values)
        IS.at[idx,"bmg"] = float((m.where(mask_shelf & (m < 0.0))*ismr2bmb*s2a).sum().values)
        
        if include_no_front:
            
            dA = dA.where(grd.mask_front == 0)
            weights = dA/dA.sum()
            
            dA_l = dA.where(m > 0.0)
            weights_l = dA_l/dA_l.sum()

            dA_g = dA.where(m < 0.0)
            weights_g = dA_g/dA_g.sum()
            
            ismr2bmb = dA*rhoi*(10**-12)
            
            IS.at[idx,"A_nf"] = float(dA.sum().values)*10**-9
            IS.at[idx,"ismr_nf"] = float((m.where(mask_shelf & (grd.mask_front == 0))*weights*s2a).sum().values)
            IS.at[idx,"ismr_nf_l"] = float((m.where(mask_shelf & (grd.mask_front == 0) & (m > 0.0))*weights_l*s2a).sum().values)
            IS.at[idx,"ismr_nf_g"] = float((m.where(mask_shelf & (grd.mask_front == 0) & (m < 0.0))*weights_g*s2a).sum().values)
            IS.at[idx,"bmb_nf"] = float((m.where(mask_shelf & (grd.mask_front == 0))*ismr2bmb*s2a).sum().values)
            IS.at[idx,"bml_nf"] = float((m.where(mask_shelf & (grd.mask_front == 0) & (m > 0.0))*ismr2bmb*s2a).sum().values)
            IS.at[idx,"bmg_nf"] = float((m.where(mask_shelf & (grd.mask_front == 0) & (m < 0.0))*ismr2bmb*s2a).sum().values)

        
    return IS%matplotlib notebook
import pandas as pd
from matplotlib.patches import Rectangle
from scipy.spatial import KDTree

def plot_ice_shelf_areas(grd,
                         google_sheet_url='https://docs.google.com/spreadsheets/d/1BpI6iRB569kF7TxdiFT3ColNVH9zwHgNwq49Dat7Ok4/edit?usp=sharing',
                         idx_start=0,idx_end=69,labels=True):

    csv_export_url = google_sheet_url.replace('edit?', 'export?gid=0&format=csv&')
    IS = pd.read_csv(csv_export_url)

    IS[['lat_min','lat_max','lon_min','lon_max']] = IS[['lat_min','lat_max','lon_min','lon_max']].convert_objects(convert_numeric=True)
    IS['eta']=np.nan
    IS['xi']=np.nan

    if labels:  
        
        lat_flat = grd.lat_rho.stack(etaxi = ('eta_rho','xi_rho'))
        lon_flat = grd.lon_rho.stack(etaxi = ('eta_rho','xi_rho'))

        points = np.column_stack((lat_flat.values,lon_flat.values))
        tree = KDTree(points)

        def find_etaxi(lat,lon,grd):

            target = np.column_stack((lat,lon))
            dist, ind = tree.query(target)

            eta = lat_flat[ind].eta_rho.values
            xi = lat_flat[ind].xi_rho.values

            return eta,xi

        for idx,row in IS.iloc[idx_start:idx_end].iterrows():
            if row.lat_min:

                IS['eta'][idx],IS['xi'][idx] = find_etaxi((row.lat_min+row.lat_max)*0.5,(row.lon_min+row.lon_max)*0.5,grd)


    plt.close()
    fig,ax = plt.subplots(figsize=(10,7))
    grd.zice.where((grd.mask_rho==1)&(grd.zice<0)).plot(ax=ax,alpha=0.2,add_colorbar=False)

    mask_ice = (grd.mask_rho==1) & (grd.zice<-10)

    for idx,row in IS.iloc[idx_start:idx_end].iterrows(): 

        if row.lon_min<row.lon_max:
            mask_coord = (grd.lon_rho>row.lon_min) & (grd.lon_rho<=row.lon_max) & (grd.lat_rho>row.lat_min) & (grd.lat_rho<=row.lat_max)
        else:
            mask_coord = ((grd.lon_rho>row.lon_min) | (grd.lon_rho<=row.lon_max)) & (grd.lat_rho>row.lat_min) & (grd.lat_rho<=row.lat_max)

        grd.zice.where(mask_ice & mask_coord).plot(ax=ax,add_colorbar=False)
        ax.contour(mask_coord,levels=1,alpha=0.5,linewidth=0.01,colors='k')
        if labels:
            ax.text(row.xi,row.eta,row.Names,fontsize=6)

    if labels:    
        latp = grd.lat_rho.plot.contour(alpha=0.2,levels=30,add_colorbar=False,colors="k")
        plt.clabel(latp,inline=1)
        lonp = grd.lon_rho.plot.contour(alpha=0.2,levels=90,add_colorbar=False,colors="k")
        plt.clabel(lonp,inline=1)

    ax.set_aspect('equal')
    plt.tight_layout()
    plt.show()plot_ice_shelf_areas(grd2)def make_mask_front(grd,nb_cells):
    
    mask_rho = grd.mask_rho.values
    mask_land = np.zeros_like(mask_rho)
    mask_land[mask_rho == 0] = 1
    mask_zice = np.zeros_like(mask_land)
    mask_zice[grd.zice.values*mask_rho != 0] = 1

    mask_front = np.zeros_like(grd.mask_rho.values)

    for j in grd.eta_rho.values:
        for i in grd. xi_rho.values:
            if mask_zice[j,i] == 1:
                j_min = max(j-nb_cells,0)
                j_max = min(j+nb_cells, np.size(mask_rho,0))
                i_min = max(i-nb_cells,0)
                i_max = min(i+nb_cells+1, np.size(mask_rho,1))

                if np.any(mask_zice[j_min:j_max,i_min:i_max] + mask_land[j_min:j_max,i_min:i_max]== 0):
                        mask_front[j,i] = 1
                        
    grd['mask_front'] = (('eta_rho','xi_rho'),mask_front)
    
    return grddef make_mass_loss_table(m,grd,
                         google_sheet_url='https://docs.google.com/spreadsheets/d/1BpI6iRB569kF7TxdiFT3ColNVH9zwHgNwq49Dat7Ok4/edit?usp=sharing',
                        idx_start=0,idx_end=69,include_no_front=True):
    
    csv_export_url = google_sheet_url.replace('edit?', 'export?gid=0&format=csv&')
    IS = pd.read_csv(csv_export_url)

    s2a = 3600*24*365
    rhoi = 916

    mask_ice = (grd.mask_rho==1) & (grd.zice<0)

    for idx,row in log_progress(IS.iloc[idx_start:idx_end].iterrows(),name='Ice shelf',every=1): 

        if row.lon_min<row.lon_max:
            mask_coord = (grd.lon_rho>row.lon_min) & (grd.lon_rho<=row.lon_max) & (grd.lat_rho>row.lat_min) & (grd.lat_rho<=row.lat_max)
        else:
            mask_coord = ((grd.lon_rho>row.lon_min) | (grd.lon_rho<=row.lon_max)) & (grd.lat_rho>row.lat_min) & (grd.lat_rho<=row.lat_max)

        mask_shelf = mask_ice & mask_coord
        
        dA = (1/(grd.pm*grd.pn)).where(mask_shelf)
        weights = dA/dA.sum()
        
        dA_l = dA.where(m > 0.0)
        weights_l = dA_l/dA_l.sum()
        
        dA_g = dA.where(m < 0.0)
        weights_g = dA_g/dA_g.sum()
        
        ismr2bmb = dA*rhoi*(10**-12)

        IS.at[idx,"A"] = float(dA.sum().values)*10**-9
        IS.at[idx,"ismr"] = float((m.where(mask_shelf)*weights*s2a).sum().values)
        IS.at[idx,"ismr_l"] = float((m.where(mask_shelf & (m > 0.0))*weights_l*s2a).sum().values)
        IS.at[idx,"ismr_g"] = float((m.where(mask_shelf & (m < 0.0))*weights_g*s2a).sum().values)
        IS.at[idx,"bmb"] = float((m.where(mask_shelf)*ismr2bmb*s2a).sum().values)
        IS.at[idx,"bml"] = float((m.where(mask_shelf & (m > 0.0))*ismr2bmb*s2a).sum().values)
        IS.at[idx,"bmg"] = float((m.where(mask_shelf & (m < 0.0))*ismr2bmb*s2a).sum().values)
        
        if include_no_front:
            
            dA = dA.where(grd.mask_front == 0)
            weights = dA/dA.sum()
            
            dA_l = dA.where(m > 0.0)
            weights_l = dA_l/dA_l.sum()

            dA_g = dA.where(m < 0.0)
            weights_g = dA_g/dA_g.sum()
            
            ismr2bmb = dA*rhoi*(10**-12)
            
            IS.at[idx,"A_nf"] = float(dA.sum().values)*10**-9
            IS.at[idx,"ismr_nf"] = float((m.where(mask_shelf & (grd.mask_front == 0))*weights*s2a).sum().values)
            IS.at[idx,"ismr_nf_l"] = float((m.where(mask_shelf & (grd.mask_front == 0) & (m > 0.0))*weights_l*s2a).sum().values)
            IS.at[idx,"ismr_nf_g"] = float((m.where(mask_shelf & (grd.mask_front == 0) & (m < 0.0))*weights_g*s2a).sum().values)
            IS.at[idx,"bmb_nf"] = float((m.where(mask_shelf & (grd.mask_front == 0))*ismr2bmb*s2a).sum().values)
            IS.at[idx,"bml_nf"] = float((m.where(mask_shelf & (grd.mask_front == 0) & (m > 0.0))*ismr2bmb*s2a).sum().values)
            IS.at[idx,"bmg_nf"] = float((m.where(mask_shelf & (grd.mask_front == 0) & (m < 0.0))*ismr2bmb*s2a).sum().values)

        
    return ISdef merge_patches(df,name1,name2,name):
    p1 = df.loc[df['Names'] == name1]
    p2 = df.loc[df['Names'] == name2]
    
    m = {'Names':name,'A':np.nan,'bmb':np.nan,'bml':np.nan,'bmg':np.nan,'ismr':np.nan,'ismr_l':np.nan,'ismr_g':np.nan,
         'A_nf':np.nan,'bmb_nf':np.nan,'bml_nf':np.nan,'bmg_nf':np.nan,'ismr_nf':np.nan,'ismr_nf_l':np.nan,'ismr_nf_g':np.nan,}
    for key in ['A','bmb','bml','bmg','A_nf','bmb_nf','bml_nf','bmg_nf']:
        m[key] = p1[key].values + p2[key].values
    for key in ['ismr','ismr_l','ismr_g','ismr_nf','ismr_nf_l','ismr_nf_g']:    
        m[key] = ((p1[key]*p1.A_nf).values + (p2[key]*p2.A_nf).values) / m['A_nf']
    
    df = pd.concat([df,pd.DataFrame(m)],ignore_index=True)

    return dfgrd2=make_mask_front(grd2,2)IS = make_mass_loss_table(avg2.m.mean('ocean_time'),grd2,include_no_front=True)ISIS_merged = merge_patches(IS,'George VI 1','George VI 2','George VI')
IS_merged = merge_patches(IS_merged,'Ronne','Filchner','Ronne-Filchner')
IS_merged = merge_patches(IS_merged,'Moscow University','Totten','Moscow Uni + Totten')
IS_merged = merge_patches(IS_merged,'Brunt/Stancomb','Riiser-Larsen','Brunt + Riiser-Larsen')
IS_merged = merge_patches(IS_merged,'Fimbul','Jelbart','Fimbul + Jelbart')
IS_merged = merge_patches(IS_merged,'Ross West','Ross East','Ross')IS_out = IS_merged[['Names','A','bmb','ismr',
                    'A_nf','bmb_nf','ismr_nf',
                    'bml','bmg','ismr_l','ismr_g',
                    'bml_nf','bmg_nf','ismr_nf_l','ismr_nf_g',
                    'lon_min','lat_min','lon_max','lat_max']] 

pd.options.display.float_format = '{:,.2f}'.format

IS_outIS_out.to_csv(csv_out_path)data = pd.read_csv(csv_out_path,index_col=0)
out_path = os.path.join(fig_dir,'mass_balance_all.png')

plt.close()
fig,ax = plt.subplots(figsize=(15,30))
data.plot(x='Names',y=['bmb','bml','bmg'],kind='barh',ax=ax,fontsize=24,colors=['green','lightgreen','lightgreen'])
ax.legend(markerscale=1,fontsize=24,loc=4)
ax.set_title('Ice shelf basal mass balance',fontsize=24)
ax.set_xlabel('Mass change rate in Gt/yr',fontsize=24)
ax.grid(axis='x')
#ax.set_xlim(-180,350)
plt.savefig(out_path,dpi=300,format='png', bbox_inches='tight')
plt.show()data = pd.read_csv(csv_out_path,index_col=0)
out_path = os.path.join(fig_dir,'ismr_balance_all.png')

plt.close()
fig,ax = plt.subplots(figsize=(15,30))
data.plot(x='Names',y=['ismr','ismr_l','ismr_g'],kind='barh',ax=ax,fontsize=24,colors=['green','lightgreen','lightgreen'])
ax.legend(markerscale=1,fontsize=24,loc=4)
ax.set_title('Ice shelf basal mass balance',fontsize=24)
ax.set_xlabel('Mass change rate in Gt/yr',fontsize=24)
ax.grid(axis='x')
#ax.set_xlim(-180,350)
plt.savefig(out_path,dpi=300,format='png', bbox_inches='tight')
plt.show()def calc_antarctic_mass_loss(m,grd):
    
    s2a = 3600*24*365.25
    rhoi = 916
    
    ice_shelf = (grd.mask_rho==1) & (grd.zice<0.0) 
    vostock = (grd.lat_rho<-75) & (grd.lat_rho>-80) & (grd.lon_rho>95) & (grd.lon_rho<115)
    
    mask = ice_shelf & ~vostock
    
    dA = (1/(grd.pm*grd.pn)).where(mask)
    weights = dA/dA.sum()
    
    ismr = (m.where(mask).mean('ocean_time')*weights*s2a).sum()
    bmb = (m.where(mask).mean('ocean_time')*dA*rhoi*s2a*10**-12).sum()
    
    print('Area of all ice shelves in 10^3 km^2: ',dA.sum().values*10**-9)
    print('Area average melt rate in m/yr: ',ismr.values)
    print('Basal mass loss in Gt/a: ',bmb.values)
    
    return%matplotlib notebook
grid_path = os.path.join(data_dir,'preprocessing','processed','waom10_grd.nc')
file_path = os.path.join(data_dir,'analysis','raw','waom10','ocean_avg_0009.nc')
grd10 = xr.open_dataset(grid_path)
avg10 = xr.open_dataset(file_path)

calc_antarctic_mass_loss(avg10.m,grd10)grid_path = os.path.join(data_dir,'preprocessing','processed','waom2_grd.nc')
file_path = os.path.join(data_dir,'analysis','raw','waom2_fix','m.nc')
grd2 = xr.open_dataset(grid_path)
avg2 = xr.open_dataset(file_path)

calc_antarctic_mass_loss(avg2.m,grd2)
# In[ ]:





# In[ ]:





# def make_mask_front(grd,nb_cells):
#     
#     mask_rho = grd.mask_rho.values
#     mask_land = np.zeros_like(mask_rho)
#     mask_land[mask_rho == 0] = 1
#     mask_zice = np.zeros_like(mask_land)
#     mask_zice[grd.zice.values*mask_rho != 0] = 1
# 
#     mask_front = np.zeros_like(grd.mask_rho.values)
# 
#     for j in grd.eta_rho.values:
#         for i in grd. xi_rho.values:
#             if mask_zice[j,i] == 1:
#                 j_min = max(j-nb_cells,0)
#                 j_max = min(j+nb_cells, np.size(mask_rho,0))
#                 i_min = max(i-nb_cells,0)
#                 i_max = min(i+nb_cells+1, np.size(mask_rho,1))
# 
#                 if np.any(mask_zice[j_min:j_max,i_min:i_max] + mask_land[j_min:j_max,i_min:i_max]== 0):
#                         mask_front[j,i] = 1
#                         
#     grd['mask_front'] = (('eta_rho','xi_rho'),mask_front)
#     
#     return grd
# grd = make_mask_front(grd,3)

# %matplotlib notebook
# s2a = 3600*24*365
# (m*s2a).where((grd.mask_rho==1) & (grd.zice<0.0)).plot(vmin=-2,vmax=2,cmap="bwr")
# grd.mask_front.plot.contour(alpha=0.3)
# plt.show()
