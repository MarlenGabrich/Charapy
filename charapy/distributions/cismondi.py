import numpy as np
from scipy.optimize import curve_fit

class Distribution_cismondi():
    def __init__(self):
        ...

    def molar_fraction(self,cn,AC,BC):
        ''' Function to calculated molar fraction

        Parameters
        ----------
        cn: int
            carbon number

        AC, BC: float
            fit parameters

        Returns
        -------
        mf: float
            fraction molar
        '''

        return np.exp(AC*cn-BC)

    def molecular_weight(self,cn,C):
        '''Functtion to calculated molecular weight

        Parameters
        ----------
        cn: int
            carbon number

        C: float
            third fit parameter - Cismondi et. al

        Returns
        -------
        mw: float
            molecular weight
        '''

        return 84+C*(cn-6)

    def density(self,cn,AD):
        '''Function to calculated density

        Parameters
        ----------
        cn: int
            carbon number
        AD: float
            fit parameter

        Returns
        -------
        rho: float
            density
        '''

        return AD*(np.exp(-cn/10)-np.exp(-0.6))+0.685

    def fit_ACBC(self,z,cn):
        '''Fit function to calculated AC, BC parameters

        Parameters
        ----------
        cn: int
            carbon number
        z: float
            molar fraction

        Returns
        -------
        AC, BC: float
            fit parameters
        '''

        dc = Distribution_cismondi()
        AB = curve_fit(dc.molar_fraction,cn,z)
        AC = AB[0][0]
        BC = AB[0][1]

        return AC,BC

    def fit_C(self,cn,mw):
        '''Fit function to calculated third parameter
        of Cismondi et. al

        Parameters
        ----------
        cn: int
            carbon number
        mw: float
            molecular weight

        Returns
        -------
        C: float
            third parameter of Cismondi et. al
        '''

        dc = Distribution_cismondi()
        C = curve_fit(dc.molecular_weight,cn,mw)[0][0]
        return C

    def fit_AD(self,cn,rho):
        '''Fit function to calculated AD parameter

        Parameters
        ----------
        cn: int
            carbon number
        rho: float
            density

        Return
        ------
        AD: float
            fit parameter
        '''

        dc = Distribution_cismondi()
        AD = curve_fit(dc.density,cn,rho)[0][0]

        return AD
