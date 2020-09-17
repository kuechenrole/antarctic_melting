import ttide as tt
import datetime
from scipy.interpolate import NearestNDInterpolator
import xarray as xr
import numpy as np
from log_progress import log_progress
import matplotlib.pyplot as plt

def NDinterp(data):

    valid_mask = ~np.isnan(data)
    coords = np.array(np.nonzero(valid_mask)).T
    values = data[valid_mask]

    it = NearestNDInterpolator(coords,values)

    filled = it(list(np.ndindex(data.shape))).reshape(data.shape)

    return filled    


def grid_ttide(da,grid_ds,stime,constit_list,res=50):
    
    ana_list = ['amp','amp_err','phase','phase_err']
    
    print('setting up the new fields ',ana_list,' for ',constit_list)
    dummy = np.empty((da.eta_rho.size,da.xi_rho.size))
    dummy[:,:] = np.nan
    
    for const in constit_list:
        for ana in ana_list:
            #print(const+'_'+ana)
            grid_ds[const+'_'+ana]=(('eta_rho','xi_rho'),dummy.copy())
     
    print("applying t_tide to every ",res,"th cell ..." )
    xi_values = np.linspace(da.xi_rho[0].values,da.xi_rho.size-1,res,dtype=int,endpoint=True)
    eta_values = np.linspace(da.eta_rho[0].values,da.eta_rho.size-1,res,dtype=int,endpoint=True)
    
    for xi in log_progress(xi_values,name='xi'):
        
        for eta in eta_values:
            da_sl = da.isel(eta_rho=eta,xi_rho=xi)
            grd_sl = grid_ds.isel(eta_rho=eta,xi_rho=xi)

            if da_sl.isnull().values.any():
                for const in constit_list:
                    for ana in ana_list:
                        grid_ds[const+'_'+ana][eta,xi]=np.NaN
                
                
            else:
                signal = da_sl.values
                latitude = grd_sl.lat_rho.values
                try:
                    ttide_out = tt.t_tide(signal,stime=stime,lat=latitude,out_style=None)
                    
                    tt_ind = {}
                    for const in constit_list:
                        tt_ind[const] = list(ttide_out['nameu']).index(str.encode(const+'  '))
                        
                        for ana,tt_ana in zip(ana_list,ttide_out['tidecon'][tt_ind[const]]):
                            grid_ds[const+'_'+ana][eta,xi] = tt_ana

                except TypeError:
                    for const in constit_list:
                        for ana in ana_list:
                            grid_ds[const+'_'+ana][eta,xi]=np.NaN
                    
    print('interpolating intermediate cells and mask land')
    for con in constit_list:
        for ana in ana_list:
            grid_ds[con+'_'+ana].values = NDinterp(grid_ds[con+'_'+ana].values)
            grid_ds[con+'_'+ana] = grid_ds[con+'_'+ana].where(grid_ds.mask_rho,0.0) 
      
        
    return grid_ds


def plot_phase(case_const_da,case_str,ref_const_da,ref_str,comp,constit):
    
    xi = []
    eta = []
    #atg_phase_diff = []
    atg_phase = []
    
    for key,sta in comp.items():
        xi.append(sta['xi_rho'])
        eta.append(sta['eta_rho'])
        atg_phase.append(sta['atg'][constit][1])
        #atg_phase_diff.append(sta['tt'][constit][2] - sta['atg'][constit][1])
    
    plt.close('all')
    fig,axes = plt.subplots(ncols=2,figsize=(15,5))
    ax1,ax2 = axes.flatten()
    
    fig.suptitle('Evaluation of '+case_str+' '+case_const_da.name+' aginst ATG and '+ref_str,fontsize=16)
    
    case_const_da.fillna(0).plot(ax=ax1,vmin=0,vmax=360)
    ref_const_da.fillna(0).plot.contour(ax=ax1,levels=np.linspace(0,360,30),linestyles='dashed',alpha=0.75)
    ax1.scatter(xi,eta,s=100,c=atg_phase,vmin=0,vmax=360,edgecolors='k')
    
    ax1.set_title('Pcolor: '+case_str+' [deg]\n  Lines: '+ref_str+' [deg]\n Scatter: ATG [deg]')
    
    phase_diff = case_const_da - ref_const_da
    pdv = phase_diff.values
    pdv[pdv>180]-=360
    pdv[pdv<-180]+=360
    phase_diff.values = pdv
    #phase_diff_rel = abs(case_const_da - ref_const_da)/360
    #atg_phase_diff_rel = np.absolute(atg_phase_diff)/360
    
    phase_diff.plot(ax=ax2)
    #ax2.scatter(xi,eta,s=100,c=atg_phase_diff_rel,vmin=0,vmax=1,edgecolors='k')
    ax2.set_title(case_str+' - '+ref_str+' difference in deg')

    for ax in axes.flatten():
        ax.set_aspect('equal')
        ax.axis("off")
        
    plt.show()
    

def plot_amp(case_const_da,case_str,ref_const_da,ref_str,comp,constit,wct_da,vmin=-0.50,vmax=0.50):
    
    xi = []
    eta = []
    atg_amp_diff = []
    
    for key,sta in comp.items():
        xi.append(sta['xi_rho'])
        eta.append(sta['eta_rho'])
        atg_amp_diff.append(sta['tt'][constit][0]-sta['atg'][constit][0])
    
    plt.close()
    fig,axes = plt.subplots(ncols=2,figsize=(15,5))
    ax1,ax2 = axes.flatten()
    
    fig.suptitle('Evaluation of '+case_str+' '+case_const_da.name+' aginst ATG and '+ref_str,fontsize=16)
    
    amp_diff = case_const_da-ref_const_da
    
    amp_diff_rel = abs(amp_diff)/ref_const_da*100
    #amp_diff_rel_norm = (abs(amp_diff)/ref_const_da*100*wct_da)/wct_da.max()
    #atg_amp_diff_rel_norm = (np.absolute(atg_amp_diff)/ref_const_da.values*100*wct_da.values)/wct_da.max().values
    
    amp_diff.plot(ax=ax1,cmap=plt.cm.bwr,vmin=vmin,vmax=vmax)
    ax1.scatter(xi,eta,s=100,c=atg_amp_diff,vmin=vmin,vmax=vmax,edgecolors='k',cmap=plt.cm.bwr)
    ax1.set_title('Pcolor: '+case_str+' - '+ref_str+' [m]\n Scatter: '+case_str+' - ATG [m]')
    
    amp_diff_rel.plot(ax=ax2,vmin=0,vmax=100)
    #ax2.scatter(xi,eta,s=100,c=atg_amp_diff_rel_norm,vmin=0,vmax=100,edgecolors='k',cmap=plt.cm.bwr)
    ax2.set_title(case_str+' - '+ref_str+' relative difference in [%]')
    
    for ax in axes.flatten():
        #ax.axis("off")
        ax.set_aspect('equal')
    
    plt.show()
