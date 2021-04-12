# equip roms grid with depths and volumes
from cartesian_grid_3d import cartesian_grid_3d
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
