import numpy as np
from scipy.optimize import curve_fit
from numpy import sqrt

class Foos():
    def __init__(self):
        self.items_list1 = []
        self.items_list2 = [] 

    def carbon_number_foo(self,z,x,y):
        """Carbon number fraction function 

        Estimates the carbon number based on a composition
        and two fit  parameters
        Use the Pedersen's ratio 
        This means that it applies from carbon 6 onward
      
        Parameters
        ----------
        z: float
          Molar fraction of component
        x, y: float
            fit parameters

        Returns 
        -------
        carbon_number: float 
            number of single carbon number 
    
        """
        for item in z:
            self.items_list1.append(item)
            self.items_list2.append[x + y*np.log(z)]

        return self.items_list2

class Foo_fit(Foos):
    def __init__(self):
        ... 

    def fit_AB(self,cn,z):
        foo = Foos.carbon_number_foo
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
        
        for item1, item2 in cn,z:
            self.items_list1.append(item1)
            self.items_list2.append(item2)

        self.A, self.B = curve_fit(foo, self.items_list1, self.items_list2)[0]

        return self.A, self.B