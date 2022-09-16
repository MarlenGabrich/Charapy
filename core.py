import numpy as np
from scipy.optimize import curve_fit
from numpy import sqrt
import pandas as pd

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
            A, B = self.A, self.B
        
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

class GeneralFit():
    def __init__(self):
        ...
    def fit_AB (self, foo, carbon_number, molar_fraction):
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
        self.A, self.B = curve_fit(foo, carbon_number,molar_fraction)[0]

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

class Correlations():
    def __init__(self,EoS):
        self.coeff = {"PR": {
            'c1':73.4043, 'c2':97.352, 'c3':0.618744, 'c4': -2059.32,
            'd1':0.0728462, 'd2':2.18811, 'd3': 163.91, 'd4': -4043.23, 'd5': 25,
            'e1': 0.373765, 'e2': 0.00549269, 'e3': 0.00117934, 'e4':-0.00000493049
            },

            "SRK": {
                'c1':163.12, 'c2':86.052, 'c3':0.43375, 'c4': -1877.4,
            'd1':-0.13408, 'd2':2.5019, 'd3': 208.46, 'd4': -3987.2, 'd5': 2.0,
            'e1': 0.74310, 'e2': 0.0048122, 'e3': 0.0096707, 'e4':-0.0000037184
            }
        }
        if isinstance(EoS, str):
            self.EoS = self.coeff[EoS]
        else:
            raise ValueError("Tiene que pasar una str")

    def coeffi(self,*kwargs):
        """ Coefficients data set to SRK or PR 

        Parameters
        ----------
        EoS: str
            Equation of state --> agregar línea para pasar el input a mayúscula para evitar errores

        Return
        ------
        c1,c2,c3,c4: float
            fit parameters to critical temperature
        d1,d2,d3,d4,d5: float
            fit parameters to critical presion
        e1,e2,e3,e4: float
            fit parameters to m_factor
        """
        x = kwargs
        if isinstance(x, str):
            return self.EoS[x]
        else:
            raise ValueError("Tiene que pasar una str válida")

    def critical_temperature(self, molecular_weight,density,
                            c1=None,c2=None,c3=None,c4=None):
        """Critical temperature correlation
        
        Parameters
        ----------
        molecular_weight: float
            Molecular weight function
        density: float
            Pedersen density
        c1,c2,c3,c4: float
            fit parameters

        Return
        ------
        critical_temperature: float
            Critical temperature
        """
        if c1 is None and c2 is None and c3 is None and c4 is None:
            c1 = self.EoS['c1']
            c2 = self.EoS['c2']
            c3 = self.EoS['c3']
            c4 = self.EoS['c4']
        
        return (
            c1*density+c2*(np.log(molecular_weight))
            +c3*molecular_weight*(c4/molecular_weight)
            )

    def critical_pression(self, molecular_weight, density, 
                        d1=None, d2=None, d3 = None, d4 = None, d5=None):
        """Critical pression correlation
        
        Parameters
        ----------
        molecular_weight: float
            Molecular weight function
        d1,d2,d3,d4,d5: float
            fit parameters

        Return
        ------
        critical_pression: float
            Critical pression
        """
        if d1 is None and d2 is None and d3 is None and d4 is None:
            d1 = self.EoS['d1']
            d2 = self.EoS['d2']
            d3 = self.EoS['d3']
            d4 = self.EoS['d4']
            d5 = self.EoS['d5']

        pression_log = d1 + d2*(density*(np.exp(d5)+(d3/molecular_weight)+(d4/(molecular_weight)*np.exp(2))))
    
        return np.exp(pression_log)
    
    def m_factor(self, molecular_weight, density, 
                e1=None, e2=None, e3=None, e4=None):

        """ *m* factor (links to acentric factor as appropiate)
        
        Parameters
        ----------
        molecular_weight: float
            Molecular weight function
        e1,e2,e3,e4: float
            fit parameters

        Return
        ------
        m_factor: float
            m factor to calculate the acentric factor
        """
        if e1 is None and e2 is None and e3 is None and e4 is None:
            e1 = self.EoS['e1']
            e2 = self.EoS['e2']
            e3 = self.EoS['e3']
            e4 = self.EoS['e4']
        
        self.m_factor = e1+e2*molecular_weight+e3*density+e4*(molecular_weight*(np.exp(2)))
        return self.m_factor

    def accentric_factor(self,m_f):
        """ Accentric factor function
        
        Parameters
        ----------
        m_factor: float
            m factor --> links to acentric factor as appropiate
        e1,e2,e3,e4: float
            fit parameters

        Return
        ------
        accentric_factor: float
            Accentric factor
        """
        if m_f is None:
            m_factor = self.m_factor

        A = 0.37464- m_factor
        B = 1.54226
        C = -0.26992

        ac_factor = (-B - sqrt(B*np.exp(2)-4*A*C))/(2*A)

        return ac_factor