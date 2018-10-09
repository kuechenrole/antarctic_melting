from netCDF4 import Dataset
from numpy import *

# Given a 2D vector in x-y space on the ROMS grid (x component on the u-grid,
# y component on the v-grid), interpolate them both to the rho-grid and rotate
# the vector to lon-lat space.
# Input:
# u = x-component of vector on the ROMS u-grid
# v = y-component of vector on the ROMS v-grid
# angle = angle between the ROMS x-axis and east, at each point, in radians
# Output:
# u_lonlat, v_lonlat = components of the vector with respect to lon-lat space
def rotate_vector_roms (u, v, angle):

    # Interpolate u to the rho-grid
    w_bdry_u = u[:,0]
    middle_u = 0.5*(u[:,0:-1] + u[:,1:])
    e_bdry_u = u[:,-1]
    u_rho = ma.concatenate((w_bdry_u[:,None], middle_u, e_bdry_u[:,None]), axis=1)
    # Interplate v to the rho-grid
    s_bdry_v = v[0,:]
    middle_v = 0.5*(v[0:-1,:] + v[1:,:])
    n_bdry_v = v[-1,:]
    v_rho = ma.concatenate((s_bdry_v[None,:], middle_v, n_bdry_v[None,:]), axis=0)

    # Rotate
    uv = (u_rho+1j*v_rho)*exp(1j*angle)
    u_lonlat = uv.real
    v_lonlat = uv.imag

    return u_lonlat, v_lonlat
    
