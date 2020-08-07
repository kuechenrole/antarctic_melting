import numpy as np
import xarray as xr
import scipy.io as sio
from scipy.spatial import KDTree
#from .log_progress import log_progress

def make_mask_sector(grd):
    mask_sector = {}
    mask_vostock = (grd.lat_rho<-75) & (grd.lat_rho>-80) & (grd.lon_rho>95) & (grd.lon_rho<115)
    mask_sector['Total Antarctica'] = xr.DataArray(np.ones_like(grd.lon_rho,dtype=bool),dims=('eta_rho','xi_rho'))
    mask_sector['Western East Antarctica'] = (grd.lon_rho>=-10.0) & (grd.lon_rho<60)
    mask_sector['Amery/Prydz Bay'] = (grd.lon_rho>=60.0) & (grd.lon_rho<80.0)
    mask_sector['Sabrina Coast/Aurora subglacial basin'] = (grd.lon_rho>=80.0) & (grd.lon_rho<130.0) & ~(mask_vostock)
    mask_sector['George V Coast/Wilkes subglacial basin'] = ((grd.lon_rho>=130.0) & (grd.lon_rho<155.0)) | ((grd.lon_rho>=155.0) & (grd.lon_rho<170.0) & (grd.lat_rho>=-72))
    mask_sector['Ross Sea'] = (grd.lon_rho>=155.0) & ~mask_sector['George V Coast/Wilkes subglacial basin'] |(grd.lon_rho<-140.0)|((grd.lon_rho>=-140.0) & (grd.lon_rho<-120.0) &(grd.lat_rho<-77.0)) 
    mask_sector['Amundsen Sea'] = ((grd.lon_rho>=-140.0) & (grd.lon_rho<-120.0) & (grd.lat_rho>=-77.0)) |((grd.lon_rho>=-120.0) & (grd.lon_rho<-90.0))
    mask_sector['Bellingshausen Sea'] = (grd.lon_rho>=-90.0) & (grd.lon_rho<-66.0) & (grd.lat_rho>=-75.0)
    mask_sector['Weddell Sea'] = (grd.lon_rho > -90) & (grd.lon_rho < -10) & ~ mask_sector['Bellingshausen Sea']
    
    return mask_sector

def make_mask_shelf_sector(grd,mask_sector,depth=1000):
    mask_shelf = ((grd.h < depth) & (grd.mask_rho ==1) & (grd.zice == 0.0)) | ((grd.zice < 0.0) & (grd.mask_rho ==1))
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
        
    return mask_shelf_sector

def make_mask_ice_shelves(coords_path,grd):
    # find nearest neighbours of ant_bounds lat_lon coords on the roms grid and 
    # define ice shelf masks and middle coordinates
    shelves = {}

    # read in ice shelf coordinates derived from ant_bounds matlab package
    coords = sio.loadmat(coords_path)['shelves'][0].squeeze()

    for idx in range(np.size(coords)):
        name = coords[idx][0][0]
        lats = coords[idx][1].squeeze()
        lons = coords[idx][2].squeeze()
        if (lats.any() == False) | (lons.any()==False):
            print('Non coordinates for '+name)
            continue
        name = coords[idx][0][0]
        shelves[name]={}
        shelves[name]['lat']= coords[idx][1].squeeze()
        shelves[name]['lon']= coords[idx][2].squeeze() 

    lat_flat = grd.lat_rho.stack(etaxi = ('eta_rho','xi_rho'))
    lon_flat = grd.lon_rho.stack(etaxi = ('eta_rho','xi_rho'))
    points = np.column_stack((lat_flat.values,lon_flat.values))
    tree = KDTree(points)

    ind_all = np.array([])
    for name,data in log_progress(shelves.items(),name='ice shelves'):

        lats = data['lat']
        lons = data['lon']

        target = np.column_stack((lats,lons))
        dist, ind = tree.query(target,k=1)

        if ind.any()==False:
            print('No grid cells found for '+name)
            continue

        ind = ind[~np.isin(ind,ind_all)]
        ind_all = np.append(ind_all,ind)

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

        shelves[name]['xi_mid']=np.median(xi)
        shelves[name]['eta_mid']=np.median(eta)

        shelves[name]['lat_mid']=np.median(lats)
        shelves[name]['lon_mid']=np.median(lons)

    return shelves                                                            
