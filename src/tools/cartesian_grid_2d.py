import numpy as np
# Given ROMS grid variables, calculates Cartesian integrands dx and dy.
# Follows the Harversine formula and linear interpolates missing
# boundary cells
# Input:
#     lon_u, lat_u, lon_v, lat_v = 2D arrays containing values for
#                                  latitude and longitude at u and
#                                  v points, with dimension
#                                  latitude(u,v) x longitude(u,v)
# Output:
#     dx, dy = 2D arrays containing Cartesian integrands with dimension
#              latitude(rho) x longitude(rho) 


def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees). Takes arrays.
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = np.absolute(lon2 - lon1) 
    dlat = np.absolute(lat2 - lat1) 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a)) 
    r = 6.371e6 # Radius of earth in meters
    
    return c * r

def cartesian_grid_2d(lon_u,lat_u,lon_v,lat_v):
    '''Calculates the cartesian distances dx, dy of computational roms grid.
    Shape(dx)=shape(dy)=shape(lon_rho)
    Usage: dx,dy = cartesian_grid_2d(lon_u,lat_u,lon_v,lat_v)'''
    
    # Calculated distances between u (530,629) and v (529,630) grid 
    dx_middle=haversine(lon_u[:,:-1],lat_u[:,:-1],lon_u[:,1:],lat_u[:,1:])
    dy_middle=haversine(lon_v[:-1,:],lat_v[:-1,:],lon_v[1:,:],lat_v[1:,:])
    
    # extrapolate missing columns or rows
    dx_west=dx_middle[:,0]+(dx_middle[:,0]-dx_middle[:,1])
    dx_east=dx_middle[:,-1]+(dx_middle[:,-1]-dx_middle[:,-2])
    
    dy_south=dy_middle[0,:]+(dy_middle[0,:]-dy_middle[1,:])
    dy_north=dy_middle[-1,:]+(dy_middle[-1,:]-dy_middle[-2,:])
    
    # merge arrays to full (530,630) shape
    dx=np.concatenate((np.expand_dims(dx_west,axis=1),dx_middle,np.expand_dims(dx_east,axis=1)),axis=1)
    dy=np.concatenate((np.expand_dims(dy_south,axis=0),dy_middle,np.expand_dims(dy_north,axis=0)),axis=0)
    
    return dx,dy
