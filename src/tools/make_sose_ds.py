import os
import sys
import xarray as xr
import pandas as pd
import scipy.io as sio
import numpy as np
import mds as mds


def make_TS_ds(sose_dir,records=None):
    '''Reads in SothernOceanStateEstimate Temperatures and Salinities and returns them in a Xarray dataset.'''
    
    #load grid data
    print("load grid")
    grid_path = os.path.join(sose_dir,'grid.mat')
    grid_raw = sio.loadmat(grid_path)
    
    longitude = grid_raw["XC"][:,0]
    latitude = grid_raw['YC'][0,:]
    depth = grid_raw['RC'].squeeze()
    DRC = grid_raw['DRC'].squeeze()
    maskC = np.flipud(np.rot90(grid_raw['maskCtrlC']))
    DXC = np.flipud(np.rot90(grid_raw['DXC']))
    DYC = np.flipud(np.rot90(grid_raw['DYC']))
    Depth = np.flipud(np.rot90(grid_raw['Depth']))
    
    #load temperature data
    print("load temperature")
    temp_path = os.path.join(sose_dir,'THETA_mnthlyBar')
    temp_raw = mds.rdmds(temp_path,100,rec=records,fill_value=np.NaN)
    
    #load salt data
    print("load salt")
    salt_path = os.path.join(sose_dir,'SALT_mnthlyBar')
    salt_raw = mds.rdmds(salt_path,100,rec=records,fill_value=np.NaN)
    
    #define array of datetime range
    time = pd.period_range('2005-01',periods=len(temp_raw),freq='M')
    time_stamp = pd.Timestamp('2005-01')
    
    #construct Xarray dataset
    print("construct Xarray dataset")
    ds = xr.Dataset({'temperature':(['time','depth','latitude','longitude'],temp_raw),
                 'salinity':(['time','depth','latitude','longitude'],salt_raw),
                 'maskC':(['latitude','longitude','depth'],maskC),
                 'DRC':(['depth'],DRC),
                 'DXC':(['latitude','longitude'],DXC),
                 'DYC':(['latitude','longitude'],DYC),
                 'Depth':(['latitude','longitude'],Depth)},
                coords={'longitude':(('longitude'),longitude),
                       'latitude':(('latitude'),latitude),
                       'depth':(('depth'),depth),
                       'time':(('time'),time),
                       'reference_time':time_stamp})
    
    print("done!")
    
    return ds
