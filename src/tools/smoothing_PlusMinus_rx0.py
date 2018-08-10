from .log_progress import log_progress
from .uvp_masks import uvp_masks
import numpy as np
import bathy_smoother

def smoothing_PlusMinus_rx0(MSK, Hobs, rx0max, AreaMatrix,rounds):
    """
    This program use the Mellor-Ezer-Oey method (Mellor et al., 1994).
    The bathymetry is optimized for a given rx0 factor by doing a sequence
    of increase/decrease at adjacent cells.

    Usage:
    RetBathy, HmodifVal, ValueFct = smoothing_PlusMinus_rx0(MSK, Hobs, rx0max, AreaMatrix)

    ---MSK(eta_rho,xi_rho) is the mask of the grid
         1 for sea
         0 for land
    ---Hobs(eta_rho,xi_rho) is the raw depth of the grid
    ---rx0max is the target rx0 roughness factor
    ---AreaMatrix(eta_rho,xi_rho) is the matrix of areas at
       rho-points.
    """

    RoughMat = bathy_smoother.bathy_tools.RoughnessMatrix(Hobs,MSK)
    print('Old Max Roughness value is: ', RoughMat.max())
    
    eta_rho, xi_rho = Hobs.shape

    ListNeigh = np.array([[1, 0],
                          [0, 1],
                          [-1, 0],
                          [0, -1]])

    RetBathy = Hobs.copy()

    HmodifVal = 0
    TheMultiplier = (1 - rx0max) / (1 + rx0max)
    tol = 0.000001
    ValueFct = 0

    for round in log_progress(range(rounds),name='steps'):
        IsFinished = 1
        for iEta in range(eta_rho):
            for iXi in range(xi_rho):
                if (MSK[iEta, iXi] == 1):
                    Area = AreaMatrix[iEta, iXi]
                    for ineigh in range(4):
                        iEtaN = iEta + ListNeigh[ineigh,0]
                        iXiN = iXi + ListNeigh[ineigh,1]
                        if (iEtaN <= eta_rho-1 and iEtaN >= 0 and iXiN <= xi_rho-1 \
                            and iXiN >= 0 and MSK[iEtaN,iXiN] == 1):
                            AreaN = AreaMatrix[iEtaN,iXiN]
                            LowerBound = RetBathy[iEtaN,iXiN] * TheMultiplier
                            if ((RetBathy[iEta,iXi] - LowerBound) < -tol):
                                IsFinished = 0
                                h = (TheMultiplier * RetBathy[iEtaN,iXiN] - RetBathy[iEta,iXi]) \
                                         / (AreaN + TheMultiplier * Area)
                                RetBathy[iEta,iXi] = RetBathy[iEta,iXi] + AreaN * h
                                RetBathy[iEtaN,iXiN] = RetBathy[iEtaN,iXiN] - Area * h
                                HmodifVal = HmodifVal + abs(h)
                                ValueFct = ValueFct + abs(h) * (Area + AreaN)
        
        if (IsFinished == 1):
            break
            
    H = AreaMatrix * Hobs * MSK
    TheBathymetry1 = H.sum()
    H = AreaMatrix * RetBathy * MSK
    TheBathymetry2 = H.sum()
    DeltaBathymetry = TheBathymetry1 - TheBathymetry2
    print('DeltaBathymetry = ', DeltaBathymetry)

    RoughMat = bathy_smoother.bathy_tools.RoughnessMatrix(RetBathy,MSK)
    print('New Max Roughness value is: ', RoughMat.max())
    
    return RetBathy, HmodifVal, ValueFct
