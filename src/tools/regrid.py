from scipy.interpolate import griddata

def regrid(sgrd,sda,tgrd):
    
#tamura comes from npstereo curvilinear grid
#we need to find points close to the 0/360 gap, wrap them and add them to the list for source points
# otherwise we get an interpolation gap between 0/360
    slon = sgrd.lon_rho.values.flatten()
    slat = sgrd.lat_rho.values.flatten()
    sdat = sda.values.flatten()
    tlon=tgrd.lon_rho.values
    tlat=tgrd.lat_rho.values

    tdat =griddata((slon,slat),sdat,(tlon,tlat),'linear')
    tda = xr.DataArray(tdat,dims=['eta_rho','xi_rho'])
    
    return tda
