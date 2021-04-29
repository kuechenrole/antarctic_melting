import os
import sys
import xarray as xr
import dask
import matplotlib.pyplot as plt
import numpy as np
import cmocean.cm as ocm
from scipy import stats
import pandas as pd
from log_progress import log_progress

def get_vrange(da,vrange):
    if vrange==None:
        mean = da.mean(dim=da.dims)
        std = da.std(dim=da.dims)
        vmin=0-2*std
        vmax=0+2*std
    elif vrange=='sat':
        vmax = max(da.max(),np.abs(da.min()))
        vmin = -vmax
    else:
        vmax=vrange[1]
        vmin=vrange[0]

    return vmin,vmax

def plot_map(da,grd,title,cbar_label,cmap,vrange=None,save=False):
    
    vmin,vmax = get_vrange(da,vrange)
    
    plt.close()
    fig,ax = plt.subplots(figsize=(10,7))
    ax.axis('off')
    ax.set_aspect('equal')
    ax.contourf(grd.mask_rho.where(grd.mask_rho==0),colors=(('0.6','0.6','0.6')))   
    plot = da.plot(ax=ax,vmin=vmin,vmax=vmax,cmap=cmap,cbar_kwargs={'label': cbar_label})
    ax.contour(grd.zice.where(grd.mask_rho==1), levels=['-0.1'], colors=('black'),linewidths=0.2)
    ax.set_title(title,fontsize=16)
    plt.tight_layout()
    if save:
        plt.savefig(os.path.join(fig_dir,'map_'+title.replace('/','_')+'.png'),format='png',bbox_inches = "tight")
    plt.show()
    
def make_is_avg(da):
    ser = pd.Series()
    for k,v in log_progress(shelves.items(),every=2):
        dA = (grd.pm*grd.pn)**-1
        weights = dA.where(v['mask'])/dA.where(v['mask']).sum()
        ser[k] = (da*weights).sum().values
    return ser.astype(float)

def make_sec_avg(da,cavity='exclude'):
    ser = pd.Series()
    for k,v in log_progress(mask_shelf.items(),every=2):
        if cavity=='only':
            mask = (v == 1) & (grd.zice < 0)
        elif cavity=='include':
            mask = (v == 1)
        elif cavity=='exclude':
            mask = (v == 1) & (grd.zice == 0)
            
        dA = ((grd.pm*grd.pn)**-1).where(mask)
        weights = dA/dA.sum()
        ser[k] = (da*weights).sum().values
    return ser.astype(float)

def make_avg_map(df,df_shelf=None,df_frame=None):
    
    v_map = np.zeros_like(grd.mask_rho)
    mask_map = np.zeros_like(v_map)
    
    for k,v in df.iteritems():
        v_map[shelves[k]['mask']]=v
        mask_map[shelves[k]['mask']]=1
    
    if isinstance(df_shelf,pd.Series):
        for k,v in df_shelf.drop('Total Antarctica').iteritems():
            mask = (mask_shelf[k] == 1) & (grd.zice == 0)
            v_map[mask]=v
            mask_map[mask]=1
            
    if isinstance(df_frame,pd.Series):
        mask_frame = np.ones_like(grd.mask_rho)
        mask_frame[50:-50,50:-50] = 0
        for k,v in df_frame.drop('Total Antarctica').iteritems():
            mask = (mask_sector[k] == 1) & (mask_frame == 1)
            v_map[mask]=v
            mask_map[mask]=1
  
    return xr.DataArray(v_map,dims=('eta_rho','xi_rho')).where(mask_map)

def make_uv_mag(u, v):

    # Interpolate u to the rho-grid
    w_bdry_u = u[:,0]
    middle_u = 0.5*(u[:,0:-1] + u[:,1:])
    e_bdry_u = u[:,-1]
    u_rho = np.ma.concatenate((w_bdry_u[:,np.newaxis], middle_u, e_bdry_u[:,np.newaxis]), axis=1)
    # Interplate v to the rho-grid
    s_bdry_v = v[0,:]
    middle_v = 0.5*(v[0:-1,:] + v[1:,:])
    n_bdry_v = v[-1,:]
    v_rho = np.ma.concatenate((s_bdry_v[np.newaxis,:], middle_v, n_bdry_v[np.newaxis,:]), axis=0)

    return xr.DataArray(xr.ufuncs.sqrt(xr.ufuncs.square(u_rho)+xr.ufuncs.square(v_rho)),dims=('eta_rho','xi_rho'))
