import numpy as np

def inverse_polar_stereo(x=None, y=None, xcentre=None, ycentre=None, ref_lat=None):
    
    #slat = -71.0
    slat = ref_lat
    slon = 0.0
    re = 6378.137
    e2 = 6.694379852e-3
    pi = 3.141592654
    e = np.sqrt(e2)

    cdr = pi / 180.0
    mflag = -1


    if (abs(slat) == 90):
        then
        
        rho = 2. * re  / np.sqrt((1 + e)  ** (1 + e)  * (1 - e)  ** (1 - e))
        
    else:
        sl = abs(slat)  * cdr
        tc = np.tan(pi / 4 - sl / 2)  / ((1 - e  * np.sin(sl))  / (1 + e  * np.sin(sl)))  ** (e  / 2)
        mc = np.cos(sl)  / np.sqrt(1 - e2  * (np.sin(sl)  ** 2))
        rho = re  * mc  / tc


    a1 = 5. / 24. * e2  ** 2
    a2 = 1. / 12. * e2  ** 3
    a3 = 7. / 48. * e2  ** 2
    a4 = 29. / 240. * e2  ** 3
    a5 = 7. / 120. * e2  ** 3


    t = np.sqrt((x - xcentre)  ** 2 + (y - ycentre)  ** 2)  / rho
    chi = (pi  / 2) - 2. * np.arctan(t)
    alat = chi + ((e2  / 2) + a1 + a2)  * np.sin(2. * chi) + (a3 + a4)  * np.sin(4. * chi) + a5  * np.sin(6. * chi)
    alat = -alat  / cdr
    alon = -np.arctan2(-(x - xcentre), (y - ycentre))  / cdr + slon
    
    return alat,alon
