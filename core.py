from tkinter import N
import numpy as np
from scipy.optimize import curve_fit

class DistributionModel:
    def __init__(self):
        """ Definite state
        ----------------------
        A,B: float
            fit parameters
        z: float
            Molar fraction of component
        """

        self.A = None
        self.B = None
        self.z = None

    def carbon_number(self,z,A,B):
        """Carbon number fraction function 

        Estimates the carbon number based on a composition
        and two fit  parameters
        Use the Pedersen's ratio 
        This means that it applies from carbon 6 onward
        
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
        if A is None and B is None:
            A,B = self.A, self.B
            return A + B * np.log(z)

    def molar_fraction(self, carbon_number, A, B): #A.self, B.self

        """ Molar fraction function

        Estimates the molar fraction based on carbon number
        and two fit  parameters
        Use the Pedersen's ratio 
        This means that it applies from carbon 6 onwards
        
        Parameters
        ----------
        carbon_number: float
          Carbon number function
        A, B: float
          fit parameters

        Returns 
        -------
        molar_fraction: float 
          Molar fraction of a given carbon number fraction 
        """
        if A is None and B is None:
            A,B = self.A, self.B
            return np.exp((carbon_number -A)/B)


class PedersenModel(DistributionModel):
    def __init__(self):
        """ Definite state
        ------------------
        L,D: float
            fit parameters
        """
        self.L = None
        self.D = None

    def density(self, carbon_number, L, D):

        """ Densities function

        Estimates the densities from C6 onwards 
        based on increase with carbon number

        Parameters
        ----------
        carbon_number: float
          Carbon number function
        L, D: float
          fit parameters

        Returns 
        -------
        density: float 
          Density of a given carbon number fraction

        """
        if L is None and D is None:
            L = self.L
            D = self.D

            return L + D*np.log(carbon_number)

    """ ¿Es necesario?
        def carbon_number(self,z):
            return self.model.carbon_number(z,self.A,self.B)
    """

    def molecular_weight(self,carbon_number):
        """ Molecular weight function

        Estimates the molecular weight based on a 
        carbon number fraction

        Parameters
        ----------
        carbon_number: float
          Carbon number function

        Returns 
        -------
        molecular_weight: float 
          Molecular weight of a given carbon number fraction 

        """
        return carbon_number*14-4 

class CismondiModel(DistributionModel):
    def __init__(self):
        """ Definite state
        ----------------------
        Ac,Bc,Ad,Bd: float
            fit parameters
        C: float
            heaviest carbon number fraction considered

        """
        Ac.self = None
        Bc.self = None
        Ad.self = None
        Bd.self = 0.685 - Ad.self*np.exp(-0.6)


    def fit(self, x, y):
        self.A, self.B = curve_fit(self.molar_fraction, x,y)[0]
        self.L, self.M = curve_fit(self.density, x, y)[0]
