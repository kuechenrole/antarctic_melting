import xarray as xr
import matplotlib.pyplot as plt
from visualization.shiftedColorMap import shiftedColorMap

def plot_ismr_map(m,mask,m_min,m_max,grd):

    m_cmap = shiftedColorMap(plt.cm.bwr,midpoint= (1 - m_max/(m_max + abs(m_min))))

    s2a = 3600*24*365.25

    plt.close()
    (m.where(mask).mean('ocean_time')*s2a).plot(vmin=m_min,vmax=m_max,size=8,cmap=m_cmap)
    plt.title('annual mean basal melt rate in [m/a]',fontsize=24)
    plt.axis('off')
    plt.show()
    
def calc_circum(m,mask,grd):
    dA = (1/(grd.pm*grd.pn)).where(mask)
    print('Area of all ice shelves in 10^3 km^2: ',dA.sum().values*10**-9)
    weights = dA/dA.sum()
    
    s2a = 3600*24*365.25

    m_avg = (m.where(mask).mean('ocean_time')*weights*s2a).sum()
    print('Area average melt rate in m/yr: ',m_avg.values)

    rhoi = 916

    bmb = m.where(mask).mean('ocean_time')*dA*rhoi*(10**-12)*s2a
    print('Basal mass loss in Gt/a: ',bmb.sum().values)
