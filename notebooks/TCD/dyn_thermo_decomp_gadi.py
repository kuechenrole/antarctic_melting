
# coding: utf-8

# In[1]:

import os
import sys
import xarray as xr
import numpy as np



proj_dir = '/home/ubuntu/bigStick/antarctic_melting'
data_dir = os.path.join(proj_dir,'data','analysis')
raw_dir = os.path.join(data_dir,'raw')
int_dir = os.path.join(data_dir,'interim')
fig_dir = os.path.join(proj_dir,'reports','tidal_melting','figures')
tab_dir = os.path.join(proj_dir,'reports','tidal_melting','tables')

src_dir = os.path.join(proj_dir,'src')
sys.path.append(src_dir)
tools_dir = os.path.join(proj_dir,'src','tools')
sys.path.append(tools_dir)


# In[12]:

file_path = os.path.join(raw_dir,'scratch','waom4','ocean_avg_hr_0010.nc')
tides = xr.open_dataset(file_path)#,chunks={'xi_rho':75,'eta_rho':75})
file_path = os.path.join(raw_dir,'scratch','waom4_nt','ocean_avg_hr_0010.nc')
no_tides = xr.open_dataset(file_path)#,chunks={'xi_rho':75,'eta_rho':75})
grid_path = os.path.join(int_dir,'grd4_dV.nc')
grd = xr.open_dataset(grid_path)

mask = xr.open_dataarray(os.path.join(int_dir,'mask_shelf.nc'))

outpath = os.path.join(int_dir,'test.nc')


# In[8]:

from dask.distributed import Client
client = Client(scheduler_file='scheduler.json')
#client = Client()
print(client)
chunks = {'xi_rho':30,'eta_rho':30}


# In[ ]:

sel = {'xi_rho' : slice(350,650),'eta_rho' : slice(700,950),
           'xi_u' : slice(350,649),'eta_u' : slice(700,950),
           'xi_v' : slice(350,650),'eta_v' : slice(700,949)}
FRIS = tides.isel(sel)
FRIS_nt = no_tides.isel(sel)
rho_sel = {'xi_rho' : slice(350,650),'eta_rho' : slice(700,950)}
FRIS_mask = ((grd.mask_rho==1) & (mask)).isel(rho_sel)

FRIS = FRIS.chunk(chunks)
FRIS_nt = FRIS_nt.chunk(chunks)

dyn = ((FRIS.ustar-FRIS_nt.ustar)*FRIS_nt.Tstar).where(FRIS_mask).mean('ocean_time')
print('ping')
therm = ((FRIS.Tstar-FRIS_nt.Tstar)*FRIS_nt.ustar).where(FRIS_mask).mean('ocean_time')
print('ping')
cov = ((FRIS.ustar-FRIS_nt.ustar)*(FRIS.Tstar-FRIS_nt.Tstar)).where(FRIS_mask).mean('ocean_time')
print('ping')
dMstar = ((FRIS.ustar*FRIS.Tstar).where(FRIS_mask).mean('ocean_time')-
          (FRIS_nt.ustar*FRIS_nt.Tstar).where(FRIS_mask).mean('ocean_time'))
print('ping')
dsum = (therm2+cov2+dyn2)


# In[15]:


xr.Dataset({'dyn':dyn,'therm':therm,'cov':cov,'dMstar':dMstar,'dsum':dsum}).to_netcdf(outpath)


# In[ ]:



