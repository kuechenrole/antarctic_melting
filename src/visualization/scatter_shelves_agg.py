import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from tools.log_progress import log_progress
import matplotlib as mpl
import cmocean.cm as ocm

mpl.rcParams.update({'font.size': 12}) 

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

def plot_regres(ax,x_all,y_all,tit,ypos,color,poly_print=False):
    xp,p,Rsquared = polyfit(x_all,y_all,2)        
    ax.plot(xp,p(xp),'-',label='constrained poly fit',alpha=0.3,c=color)
    if poly_print:
        ax.text(0.01,ypos,tit+':'+str(p)+' ($r^2$ = %.2f)' %(Rsquared),transform=ax.transAxes,fontsize=12,color=color)
    else:
        ax.text(0.01,ypos,tit+': ($r^2$ = %.2f)' %(Rsquared),transform=ax.transAxes,fontsize=12,color=color)
        
    return p
    
def scatter_shelves_agg(x_quant,y_quant,title,xlab,ylab,shelves_dict,save=False):
    
    n = len(shelves_dict)
    colors = iter(ocm.phase(np.linspace(0,1,n)))
    shelfcol = {} 
    for key,data in shelves_dict.items():
        shelves_dict[key]['color'] = next(colors)
    
    plt.close()
    fig,ax = plt.subplots(figsize=(12,10))
    
    x_all=[]
    y_all=[]

    ypos=0.95
    for key,data in log_progress(shelves_dict.items(),every=2):
        
        if data['A']<5: continue
        
        c=data['color']
            
        x = x_quant.squeeze().where(data['mask']).values.flatten()
        y = y_quant.squeeze().where(data['mask']).values.flatten()
        mask = (np.isnan(x)==False)
        mask[np.isnan(y)]=False
        x = x[mask]
        y = y[mask]
        
        #ax.plot(x,y,'.',ms=0.1,c=c,label=key)
        
        xp,p,Rsquared = polyfit(x,y,2)
        ax.plot(xp,p(xp),c=c,label='linear fit',alpha=0.3)
            
        shelves_dict[key]['polyfit'] = p
        shelves_dict[key]['r2'] = Rsquared
        
        plot_regres(ax,x,y,key,ypos,c)
        ypos=ypos-0.025
       # x_all.extend(x)
       # y_all.extend(y)
        
    for key,data in log_progress(shelves_dict.items(),every=2):
        
        if data['A']<5:
            c='k'
            ms='.'
        else:
            c=data['color']
            ms='o'
        
        x = x_quant.squeeze().where(data['mask']).mean().values
        y = y_quant.squeeze().where(data['mask']).mean().values
        
        x_all.append(x)
        y_all.append(y)
        
        ax.plot(x,y,ms,markeredgecolor='k',label=key,c=c)
    
        
        
        if data['A']>5:
            #lf = "%.2f" %(data['r2'])
            ax.annotate(key,(x,y),size='xx-small')
        
    x_all = np.array(x_all).squeeze()
    y_all = np.array(y_all).squeeze()
    
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    ax.set_xlim([min(x_all)-0.5*np.std(x_all),max(x_all)+0.5*np.std(x_all)])
    ax.set_ylim([min(y_all)-0.5*np.std(y_all),max(y_all)+0.5*np.std(y_all)])
    
    plot_regres(ax,x_all,y_all,'ice shelf averages',0.975,'k')
    
    if save:
        plt.savefig(os.path.join(fig_dir,'scatter_shelves_agg_'+title.replace('/','_')+'.png'),format='png',bbox_inches = "tight")    

    plt.show()   