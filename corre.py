import numpy as np
from numpy import sqrt

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