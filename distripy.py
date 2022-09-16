import numpy as np
from fity import GeneralFit 

class DistributionModel():
    def __init__(self):
        ...
        
    def carbon_number(self, z, A=None, B=None):
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
        if A or B is None:
            A,B = self.A, self.B
        
        return (A+B*np.log(z))

class PedersenModel(DistributionModel):
    def __init__(self):
        ...

    def molar_fraction(self, carbon_number, A=None, B=None):
        """Molar fraction function
        
        Estimates the molar fraction based on carbon number
        and two fit parameters
        Use the Pedersen's ratio
        This means that it applies from carbon 6 onwards
        
        Parameters
        ----------
        carbon_number: float
            Carbon number function
        A,B: float
            fit parameters

        Returns
        -------
        molar_fraction: float
            Molar fraction of a given carbon number fraction
        """
        if A is None and B is None:
            A,B = self.A, self.B
    
        return np.exp(A+B*carbon_number)

    def molecular_weight(self, carbon_number):
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
    
    def density(self, carbon_number, L=None, D=None):
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

class CismondiModel(DistributionModel):
    def __init__(self):
        ...
    
    def molar_fraction(self, carbon_number, Ac=None, Bc=None):
        """ Molar fraction function

        Estimates the molar fraction based on carbon number
        and two fit  parameters
        
        Parameters
        ----------
        carbon_number: float
          Carbon number function
        Ac, Bc: float
          fit parameters

        Returns 
        -------
        molar_fraction: float 
          Molar fraction of a given carbon number fraction 
        """

        if Ac is None and Bc is None:
            Ac,Bc = self.Ac, self.Bc
      
        return np.exp(Ac*carbon_number + Bc)

    def molecular_weight(self,carbon_number,C=None):
    
        """ Molecular weight function

        Estimates the molecular weight based on carbon number
        and the third fit parameter of Cismondi  "C"

        Parameters
        ----------
        carbon_number: float
            Carbon number function
        c: float
            fit parameter (third parameter of cismondi et al.)

        Returns 
        -------
        molecular_weight: float 
            Molecular weight of a given carbon number fraction
        """

        if C is None:
            C = self.C

        return 84 + C*(carbon_number-6)

    def density(self,carbon_number,Ad=None):
        """ Density function
        Estimates the densities from C6 onwards 
        based on [...]
      
        Parameters 
        ----------
        carbon_number: float
            Carbon number function
        Ad: float
            fit parameters of cismondi's propose

        Returns
        -------
        density: float
            Density of a given carbon number fraction
        """

        if Ad is None:
            Ad = self.Ad

        Bd = 0.685 - Ad*np.exp(-0.6)
        return Ad*(np.exp(-carbon_number/10))+Bd