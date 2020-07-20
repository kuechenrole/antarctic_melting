import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os

proj_dir = os.path.join(os.pardir)
data_dir = os.path.join(proj_dir,'data')

grd_path = os.path.join(data_dir,'preprocessing','processed','waom10_grd.nc')
grd = xr.open_dataset(grd_path)

mask = (grd.zice<0) & (grd.mask_rho==1)

his_path = os.path.join(data_dir,'analysis','raw','scratch','waom10','ocean_his.nc')
m_his = xr.open_dataset(his_path).m.where(mask).mean({'xi_rho','eta_rho'})

avg_path = os.path.join(data_dir,'analysis','raw','scratch','waom10','ocean_avg.nc')
m_avg = xr.open_dataset(avg_path).m.where(mask).mean({'xi_rho','eta_rho'})


m_avg.plot(label='avg')
m_his.plot(label='his')
plt.legend()
plt.savefig('his_avg_test.png')
