import numpy as np
import xarray as xr
from .log_progress import log_progress

def dTemp(Sal,Temp,Pres):
    
    '''Calculates from the salinity (Sal,psu), the in-situ Temperature
       (Temp, degC) and the in-situ pressure (Pres, dbar) the adiabatic 
       temperature gradient (dTemp, K Dbar^-1).

         Check values:  dTemp  =     3.255976E-4 K dbar^-1
                  given Sal    =    40.0         psu
                        Temp   =    40.0         degC
                        Pres   = 10000.000       dbar
    ''' 
    s0=35.0
    a0,a1,a2,a3 = (3.5803e-5, 8.5258e-6, -6.8360e-8, 6.6228e-10)
    b0,b1 = (1.8932e-6, -4.2393e-8)
    c0,c1,c2,c3 = (1.8741e-8, -6.7795e-10, 8.7330e-12, -5.4481e-14)
    d0,d1  = (-1.1351e-10, 2.7759e-12)
    e0,e1,e2 = (-4.6206e-13,  1.8676e-14, -2.1687e-16)
    ds = Sal-s0

    dTemp = (((e2*Temp + e1)*Temp + e0 )*Pres + \
            ((d1*Temp + d0)*ds + ((c3*Temp + c2)*Temp + c1 )*Temp + c0 ) )*Pres + \
            (b1*Temp + b0)*ds +  ((a3*Temp + a2)*Temp + a1 )*Temp + a0

    return dTemp


def thetaa(Sal,Temp,Pres,RPres):
    '''! Calculates from the salinity (sal, psu), the in-situ temperature 
    ! (Temp, degC) and the in-situ pressure press, dbar) the potential 
    ! temperature (Theta, degC) converted to the reference pressure
    ! (RPres, dbar). A Runge-Kutta procedure of the fourth order is used.
    !
    ! Check value: theta   =    36.89073  degC
    !         given sal    =    40.0      psu
    !               Temp   =    40.0      degC
    !               pres   = 10000.000    dbar
    !               rfpres =     0.000    dbar'''

    ct2,ct3 =  (0.29289322 ,  1.707106781)
    cq2a,cq2b = (0.58578644 ,  0.121320344)
    cq3a,cq3b = (3.414213562, -4.121320344)
    p  = Pres
    t  = Temp
    dp = RPres-Pres
    dt = dp*dTemp(Sal,t,p)
    t  = t +0.5*dt
    q = dt
    p  = p +0.5*dp
    dt = dp*dTemp(Sal,t,p)
    t  = t + ct2*(dt-q)
    q  = cq2a*dt + cq2b*q
    dt = dp*dTemp(Sal,t,p)
    t  = t + ct3*(dt-q)
    q  = cq3a*dt + cq3b*q
    p  = RPres
    dt = dp*dTemp(Sal,t,p)
    thetaa = t + (dt-q-q)/6.0

    return thetaa

def potit(Sal,theta,Pres,RPres):
    '''! *********************************************************************
    ! Calculates from the salinity (sal, psu), potential temperature 
    ! (theta, degC) and reference pressure (pres, dbar) the in-situ 
    ! temperaure (Temp_insitu, degC) related to the in-situ pressure 
    ! (rfpres, dbar) with the help of an iterative method.'''

    tpmd = 0.001
    epsi = 0.0
    for ind in np.arange(1,100):
        Temp   = theta+epsi
        thetad  = thetaa(Sal,Temp,Pres,RPres)-theta
        if(np.abs(thetad) < tpmd):
            return Temp
        epsi = epsi-thetad

    print(' WARNING! in-situ temperature calculation has not converged!')

def calc_insituTemp(mask,temp,salt,zice):
    situTemp = temp.where(mask).copy()
    situTemp.attrs['long_name']= 'time-averaged insitu temperature at ice base'
    dummy = np.full(situTemp.values.shape,np.nan)

    for eta in log_progress(temp.eta_rho):
        if mask.sel(eta_rho=eta).any():
            for xi in temp.xi_rho:
                if mask.sel(xi_rho=xi,eta_rho=eta):
                    sali = salt.sel(xi_rho=xi,eta_rho=eta).values
                    tempi = temp.sel(xi_rho=xi,eta_rho=eta).values
                    presi = -zice.sel(xi_rho=xi,eta_rho=eta).values
                    situTempi = potit(sali,tempi,presi,0)
                    #print(sali,tempi,presi,situTempi)
                    dummy[eta,xi] = situTempi

    situTemp.values = dummy
    return situTemp

def calc_ustar(u,v,mask):
    
    Cd=5.0e-3
    Cdrt = np.sqrt(Cd)
    u = u.isel(s_rho=30).values
    v = v.isel(s_rho=30).values
    # Interpolate u to the rho-grid
    w_bdry_u = u[:,0]
    middle_u = 0.5*(u[:,0:-1] + u[:,1:])
    e_bdry_u = u[:,-1]
    u_rho = np.concatenate((w_bdry_u[:,None], middle_u, e_bdry_u[:,None]), axis=1)
    # Interplate v to the rho-grid
    s_bdry_v = v[0,:]
    middle_v = 0.5*(v[0:-1,:] + v[1:,:])
    n_bdry_v = v[-1,:]
    v_rho = np.concatenate((s_bdry_v[None,:], middle_v, n_bdry_v[None,:]), axis=0)
    
    mag = np.sqrt(u_rho**2+v_rho**2)
    ustar = Cdrt*mag
    return xr.DataArray(ustar,dims=('eta_rho','xi_rho'),name='ustar').where(mask)

def calc_frc(Tm,Sm,zice,ustar,f,m,mask):
    
    Pr = 13.8
    Sc = 2432.2
    L = 3.33e5
    cp_w = 3947.0
    rho_i = 920.0
    Ti = -20.0
    Si = 0.0
    Pradj = 12.5*(Pr**(2.0/3.0))
    cp_i =  152.5+7.122*(273.15+Ti)
    a = -0.057
    b = 0.0939
    c = 7.61e-4
    
    TFb = a*Sm + b + c*zice
    Tstar = Tm - TFb
    mflag = xr.ufuncs.sign(Tm-TFb)
    mflag=xr.ufuncs.sign((mflag+1)*mflag)
    TFi = (1 - mflag)*a*Si + b + c*zice
    turb = 2.5*np.log(5300.0*ustar*ustar/np.abs(f))
    gammaT = ustar/(turb + Pradj - 6.0)
    Tb = (gammaT*Tm+mflag*(cp_i/cp_w)*m*Ti-(L/cp_w)*m)/(gammaT + mflag*(cp_i/cp_w)*m)
    Sb = (Tb - b - c*zice)/a
    
    Tdr = xr.DataArray(Tm - TFb,dims=('eta_rho','xi_rho'),name='Tdr').where(mask)
    Tfr = xr.DataArray(Tm-Tb,dims=('eta_rho','xi_rho'),name='Tfr').where(mask)
    
    return Tdr,Tfr
