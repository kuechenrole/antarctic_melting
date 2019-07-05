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

def plot_regres(ax,x_all,y_all,tit,ypos,color,poly_print=False):
    xp,p,Rsquared = polyfit(x_all,y_all,2)        
    ax.plot(xp,p(xp),'-',label='constrained poly fit',alpha=0.3,c=color)
    if poly_print:
        ax.text(0.01,ypos,tit+':'+str(p)+' ($r^2$ = %.2f)' %(Rsquared),transform=ax.transAxes,fontsize=12,color=color)
    else:
        ax.text(0.01,ypos,tit+': ($r^2$ = %.2f)' %(Rsquared),transform=ax.transAxes,fontsize=12,color=color)
        
    return p

def scatter_shelves_avg(x_name,y_name,title,xlab,ylab,shelves_dict=shelves2,sector_dict=sector2,big=False,save=False):
    
    matplotlib.rcParams.update({'font.size': 18})
    
    plt.close()
    fig,ax = plt.subplots(figsize=(12,10))
    
    ypos=0.87
    
    x_all = []
    y_all = []
    
    for sec_key,sec_data in sector_dict.items():
        
        if sec_key=='Total Antarctica': continue
        
        x_sec = []
        y_sec = []
        
        shelves_dict_sel = {k:v for k,v in shelves_dict.items() if v['sector']==sec_key}

        for shelf_key,shelf_data in shelves_dict_sel.items():

            if big:
                if shelf_data['A'] > 5:
                    continue

            x = shelf_data[x_name]
            y = shelf_data[y_name]

            ax.plot(x,y,'.',c=sec_data.color)
            #        ax.annotate(key,(x,y),size='x-small')

            x_sec.append(x)
            y_sec.append(y)
            
            x_all.append(x)
            y_all.append(y)
    
        x_sec = np.array(x_sec).squeeze()
        y_sec = np.array(y_sec).squeeze()
    
        plot_regres(ax,x_sec,y_sec,sec_key,ypos,sec_data.color,poly_print=True)    
        ypos=ypos-0.06
    
    x_all = np.array(x_all).squeeze()
    y_all = np.array(y_all).squeeze()
    
    plot_regres(ax,x_all,y_all,'Total Antarctica',0.94,'k',poly_print=True)
    
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title,fontsize=18)
    
    if save:
        plt.savefig(os.path.join(fig_dir,'tmp','scatter_shelfavg_'+title.replace('/','_').replace(' ','_')+'.png'),format='png',dpi=300,bbox_inches = "tight")    
    
    plt.show()
    
    return