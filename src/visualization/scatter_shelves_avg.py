import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from tools.log_progress import log_progress

# Polynomial Regression
def polyfit(x, y, degree, zeroCept=False):
    results = {}

    if zeroCept==False:
        coeffs = np.polyfit(x, y, degree)
    else:
        coeffs = np.append(np.polyfit(x, y/x, degree-1),0)

    # r-squared
    p = np.poly1d(coeffs)
    # fit values, and mean
    yhat = p(x)                         # or [p(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    Rsquared = ssreg / sstot
    
    xp =np.linspace(min(x),max(x),100)
    
    return xp,p,Rsquared

def find_sector(shelves_dict,sector_dict):
    
    n = len(sector_dict)-1
    colors = iter(plt.cm.rainbow(np.linspace(0,1,n)))
    seccol = {} 
    for key,data in sector_dict.items():
        if key == 'Total Antarctica': continue
        sector_dict[key].attrs['color'] = next(colors)
    
    for key_shelf,data_shelf in shelves_dict.items():
        mask_shelf = data_shelf['mask']
        nb_shelf = np.count_nonzero(mask_shelf)
        for key_sec,mask_sec in sector_dict.items():
            if key_sec == 'Total Antarctica': continue
            nb_overlap = np.count_nonzero(mask_sec & mask_shelf)
            if nb_shelf-nb_overlap < (nb_shelf/2):
                shelves_dict[key_shelf]['sector']=key_sec
                shelves_dict[key_shelf]['sector_color']=sector_dict[key_sec].attrs['color']
                continue
                
       
    return shelves_dict,sector_dict


def make_shelves_avg(da,shelves_dict,name):
    for k,v in log_progress(shelves_dict.items(),every=2):       
        if v['mask'].any()==False:
            print('No mask for '+k)
            continue
        shelves_dict[k][name] = da.squeeze().where(v['mask']).mean().values
    return shelves_dict