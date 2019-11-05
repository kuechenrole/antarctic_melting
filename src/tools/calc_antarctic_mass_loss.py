
import xarray as xr

def calc_antarctic_mass_loss(m,grd):
    
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

    return dA.sum().values*10**-9,ismr.values,bmb.values
    
