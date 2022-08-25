from tkinter import N
import numpy as np
from scipy.optimize import curve_fit

class DistributionModel:
    def __init__(self):
        ...
    
    def carbon_number(self,z,A,B):
        """Carbon number function

        Estimates the carbon number based on a composition
        and two fit  parameters

        Parameters
        ----------
        z: float
          Molar fraction of component
        A, B: float
          fit parameters

        Returns 
        -------
        carbon_number: float 
          number of single carbon number 
    
        """
    def molar_fraction(self, carbon_number, A=None, B=None):
        if A is None and B is None:
            A,B = self.A, self.B
            return np.exp((carbon_number -A)/B)

class PedersenModel(DistributionModel):
    def __init__(self):
        self.A = None
        self.B = None

        self.L = None
        self.M = None

    def density(self, carbon_number, L=None, M=None):
        if L is None and M is None:
            L = self.L
            M = self.M

            return L + M*np.log(carbon_number)

    def carbon_number(self,z):
        return self.model.carbon_number(z,self.A,self.B)
    
    def fit(self, x, y):
        self.A, self.B = curve_fit(self.molar_fraction, x,y)[0]
        self.L, self.M = curve_fit(self.density, x, y)[0]
