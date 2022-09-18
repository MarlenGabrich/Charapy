#!usr/bin/env python3
from typing import List
import numpy as np
from math import sqrt
from scipy.optimize import curve_fit

class Foos:
    def __init__(self):
        print('in init')
        ...

    def carbon_number_foo(self,z,*A,B) -> dict or list or float:
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
        return A + B*np.log(z)

class Distribution_pedersen (Foos):
    def __init__(self):
        ...
    def p_density(self, cn, *L, D):
        """ Densities function

        Estimates the densities from C6 onwards 
        based on increase with carbon number

        Parameters
        ----------
        cn: float
            Carbon number function
        L, D: float
            fit parameters

        Returns 
        -------
        density: float 
            Density of a given carbon number fraction

        """
        return L + D*np.log(cn)
    def p_molar_fraction(self, cn, *A, B):
        """Molar fraction function
        
        Estimates the molar fraction based on carbon number
        and two fit parameters
        Use the Pedersen's ratio
        This means that it applies from carbon 6 onwards
        
        Parameters
        ----------
        cn: float
            Carbon number function
        A,B: float
            fit parameters

        Returns
        -------
        molar_fraction: float
            Molar fraction of a given carbon number fraction
        """
    
        return np.exp(A + B*cn)
    def p_molecular_weight(self, cn):
        """ Molecular weight function

        Estimates the molecular weight based on a 
        carbon number fraction

        Parameters
        ----------
        cn: float
            Carbon number function

        Returns 
        -------
        molecular_weight: float 
            Molecular weight of a given carbon number fraction 

        """
        return cn*14-4

class Distribution_cismondi(Foos):
    def __init__(self):
        ...
    def c_molar_fraction(self, cn, *Ac, Bc):
        """ Molar fraction function

        Estimates the molar fraction based on carbon number
        and two fit  parameters
        
        Parameters
        ----------
        cn: float
          Carbon number function
        Ac, Bc: float
          fit parameters

        Returns 
        -------
        molar_fraction: float 
          Molar fraction of a given carbon number fraction 
        """

        return np.exp(Ac*cn + Bc)

    def c_molecular_weight(self,cn,*C):
    
        """ Molecular weight function

        Estimates the molecular weight based on carbon number
        and the third fit parameter of Cismondi  "C"

        Parameters
        ----------
        cn: float
            Carbon number function
        C: float
            fit parameter (third parameter of cismondi et al.)

        Returns 
        -------
        molecular_weight: float 
            Molecular weight of a given carbon number fraction
        """

        return 84 + C*(cn-6)

    def c_density(self,cn,*Ad):
        """ Density function
        Estimates the densities from C6 onwards 
        based on [...]
      
        Parameters 
        ----------
        cn: float
            Carbon number function
        Ad: float
            fit parameters of cismondi's propose

        Returns
        -------
        density: float
            Density of a given carbon number fraction
        """

        Bd = 0.685 - Ad*np.exp(-0.6)
        return Ad*(np.exp(-cn/10))+Bd

class Foo_fit(Distribution_pedersen, Distribution_cismondi): 
    def __init__(self):
        self.foos = Foos()

    def fit_AB(self,cn,z):
        """ Parameter setting function

        Parameters
        ----------
        cn: float
            carbon number
        z: float
            molar fraction

        Return
        ------
        A,B: float
            fit parameters
        """

        self.A, self.B = curve_fit(self.foos.carbon_number_foo, cn, z, p0=[1,1])[0]

        return {'A':self.A, 'B': self.B}

    def fit_LM(self, cn, d):
        """ Parameter setting function

        Parameters
        ----------
        cn: float
            carbon number
        d: float
            density

        Return
        ------
        L,M: float
            fit parameters
        """ 
        self.L, self.M = curve_fit(Distribution_pedersen().p_density, cn, d)[0]

        return {'L': self.L, 'M': self.M}
    
    def fit_C(self,cn,mw):
        """ Parameter setting function
        Parameters
        ----------
        x: float
            carbon number 
        mw: float
            molecular_weight
        Returns
        -------
        C: float
            Third parameter of Cismondi's model
        """ 
        self.C = curve_fit(Distribution_cismondi().c_molecular_weight,cn,mw)[0]

        return self.C
    
    def fit_Ad(self,x,d):
        """ Parameter setting function
        Parameters
        ----------
        x: float
            carbon number 
        d: float
            density for Ad
        Returns
        -------
        C: float
            Third parameter of Cismondi's model
        """ 
        
        self.Ad = curve_fit(Distribution_cismondi().c_density,x,d)[0]
        return {'Ad': self.Ad}

    def fit_AcBc(self,x,mf):
        """ Parameter setting function
        Parameters
        ----------
        x: float
            carbon number 
        d: float
            molar_fraction for Ac, Bc
        Returns
        -------
        Ac, Ad: float
            fit parameter
        """ 
        self.Ac, self.Bc = curve_fit(Distribution_cismondi().c_molar_fraction,x,mf)[0]

        return {'Ac': self.Ac, 'Bc': self.Bc}   

class Correlations:
    def __init__(self,EoS,*x):
        self.dic_coeff = {
        "PR": {'c1':73.4043, 'c2':97.352, 'c3':0.618744, 'c4': -2059.32,
            'd1':0.0728462, 'd2':2.18811, 'd3': 163.91, 'd4': -4043.23, 'd5': 25,
            'e1': 0.373765, 'e2': 0.00549269, 'e3': 0.00117934, 'e4':-0.00000493049},

        "SRK": {
                'c1':163.12, 'c2':86.052, 'c3':0.43375, 'c4': -1877.4,
            'd1':-0.13408, 'd2':2.5019, 'd3': 208.46, 'd4': -3987.2, 'd5': 2.0,
            'e1': 0.74310, 'e2': 0.0048122, 'e3': 0.0096707, 'e4':-0.0000037184}
            }

        if isinstance(EoS, str):
            if x and isinstance(x, str): 
                return self.dic_coeff[EoS](x)
            return self.dic_coeff[EoS]
        else:
            raise ValueError("Tiene que pasar una str")
            
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