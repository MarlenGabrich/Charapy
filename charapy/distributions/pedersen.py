from scipy.optimize import curve_fit
import numpy as np

class Distribution_pedersen():
    '''Estimation ratio of Pedersen can used about
    carbon number 6 onwards'''

    def __init__(self):
        ...
    
    def carbon_number(self,z,A,B):
        '''General function to estimated carbon number
        based on components fraction molar

        Parameters
        ----------
        z: float
            component fraction molar
        A,B: float 
            fit parameters
        
        Returns
        -------
        cn: int
            carbon number
        '''

        cn = A*np.log(z) + B
        return cn

    def density(self,cn,L,M):
        '''Function to calculated density based on
        carbon number
        
        Parameters
        ----------
        cn: int
            carbon number
        L,D: float
            fit parameters

        Returns
        -------
        rho: float
            density of given carbon number
        '''

        rho = L*np.log(cn) + M
        return rho

    def molar_fraction(self,cn,A,B):
        '''Funtion to calculated molar fraction

        Parameters
        ----------
        cn: float
            carbon_number
        A,B: float
            fit parameters

        Returns
        -------
        z: float
            molar fraction
        '''

        return np.exp((cn-B)/A)

    def molecular_weight(self,cn):
        '''Function to calculated molecular weight

        Parameters
        ----------
        cn: int
            carbon number

        Returns
        -------
        mw: float
            molecular weight
        '''

        return cn*14-4

    def fit_AB(self,z,cn):
        '''Fit function to A, B parameters

        Parameters
        ----------
        cn: int
            carbon number
        z: float
            molar fraction

        Return
        ------
        A,B: float
            fit parameters and variance
            covariance matrix
        '''
        dp = Distribution_pedersen()
        AB = curve_fit(dp.carbon_number,z,cn)
        A = AB[0][0]
        B = AB[0][1]

        return A,B

    def fit_LM(self,cn,density):
        '''Fit function to L,M parameters
        
        Parameters
        ----------
        cn: int
            carbon number
        density: float
            density

        Returns
        -------
        L,M: float
            fit parameters and variance
            covariance matrix
        '''

        dp = Distribution_pedersen()
        LM = curve_fit(dp.density,cn,density)
        L = LM[0][0]
        M = LM[0][1]

        return L,M