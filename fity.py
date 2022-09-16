from distripy import DistributionModel, PedersenModel, CismondiModel
import numpy as np
from scipy.optimize import curve_fit

class GeneralFit(DistributionModel):
    def __init__(self):
        ...
    def fit_AB (self, carbon_number, molar_fraction):
        """ Parameter setting function

        Parameters
        ----------
        carbon_numer: float
            carbon number
        molar_fraction: float
            molar fraction

        Return
        ------
        A,B: float
            fit parameters
        """
        self.A, self.B = curve_fit(self.carbon_number, carbon_number,molar_fraction)[0]

        return self.A, self.B
        
class Pedersen_Fit(PedersenModel, GeneralFit):
    def __init__(self):
        ...
    def fit(self,x, y = None):
        """ Parameter setting function

        Parameters
        ----------
        x: float
            carbon number
        y: float
            molar_fraction

        Return
        ------
        L,M: float
            fit parameters
        """ 
        self.L, self.M = curve_fit(self.density, x)[0]

        if y:
            A,B = self.A, self.B

        return {'A': A,'B': B,'L': self.L, 'M': self.M}

class Cismondi_Fit(CismondiModel, GeneralFit):
    def __init__(self):
        ...
    def fit(self,x, l = None, m = None, p = None):

        """ Parameter setting function

        Parameters
        ----------
        x: float
            carbon number 

        l: float
            molecular_weight for C
        m: float
            density for Ad
        p: float
            molar_fraction for Ac, Bc
        
        Returns
        -------
        C: float
            Third parameter of Cismondi's model
        Ac, Ad: float
            fit parameter
        """ 
        self.C = '-'
        self.Ad = '-'
        self.Ac, self.Bc = '-'

        if l:
            self.C = curve_fit(CismondiModel.molecular_weight,x,l)[0]
        if m:
            self.Ad = curve_fit(CismondiModel.density,x,m)[0]
        if p:
            self.Ac, self.Bc = curve_fit(CismondiModel.molar_fraction,x,p)[0]

        return {'C': self.C, 'Ac': self.Ac, 'Ad': self.Ad}