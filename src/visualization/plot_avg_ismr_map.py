import xarray as xr
import cmocean.cm as ocm
import matplotlib.pyplot as plt
from visualization.shiftedColorMap import shiftedColorMap
from matplotlib import rcParams

def plot_avg_ismr_map(m,grd,m_min=-1,m_max=8,out_path=None):


    m_cmap = shiftedColorMap(ocm.balance,midpoint= (1 - m_max/(m_max + abs(m_min))))

    s2a = 3600*24*365.25
    mask = (grd.zice < 0.0)&(grd.mask_rho==1)
    
    land_zice = (grd.mask_rho==0) | (grd.zice < 0.0)
    plt.close()
    fig,ax = plt.subplots(figsize=(10,8))
    ax.contourf(grd.mask_rho.where(land_zice).values,colors=(('0.6','0.6','0.6')))
    (m.where(mask).mean('ocean_time')*s2a).plot(vmin=m_min,vmax=m_max,ax=ax,cmap=m_cmap)
    ax.contour(-grd.zice.where(grd.mask_rho).values, levels=['0.01'], colors=('black'),linewidths=0.5)
    ax.contour(grd.mask_rho.values, levels=['0.01'], colors=('black'),linewidths=0.5)
    plt.title('Ice shelf melt rate (m/y), annual average', fontsize=20)
    ax.set_aspect('equal')
    ax.axis('off')
    if out_path:
        plt.savefig(out_path)
    plt.show()
    
    

def plot_avg_ismr_diff_map(m_from,m_to,grd,m_min=-1,m_max=1,p_min=-1000,p_max=1000,out_path=None):

    rcParams.update({'font.size': 18})
    s2a = 3600*24*365.25
    mask = (grd.zice < 0.0)&(grd.mask_rho==1)
    land_zice = (grd.mask_rho==0) | (grd.zice < 0.0)

    m_diff = (m_to-m_from)*s2a
    m_diff_rel = (m_to-m_from)/m_from*100

    plt.close()
    fig,axes = plt.subplots(ncols=2,figsize=(20,8))
    ax1,ax2 = axes.flatten()
    for ax in axes:
        ax.contourf(grd.mask_rho.where(land_zice).values,colors=(('0.6','0.6','0.6')))
        ax.contour(-grd.zice.where(grd.mask_rho).values, levels=['0.01'], colors=('black'),linewidths=0.5)
        ax.contour(grd.mask_rho.values, levels=['0.01'], colors=('black'),linewidths=0.5)
        ax.set_aspect('equal')
        ax.axis('off')


    (m_diff.where(mask).mean('ocean_time')).plot(vmin=m_min,vmax=m_max,ax=ax1,cmap=plt.cm.seismic)
    ax1.set_title('Melt rate difference (m/y), annual average', fontsize=20)
    (m_diff_rel.where(mask).mean('ocean_time')).plot(vmin=p_min,vmax=p_max,ax=ax2,cmap=plt.cm.seismic)
    ax2.set_title('Relative melt rate difference (%), annual average', fontsize=20)
    plt.tight_layout()
    if out_path:
        plt.savefig(out_path,bbox_tight=True,dpi=300)

    plt.show()
