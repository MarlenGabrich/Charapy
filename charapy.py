import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


class Foos():
    com = {'cn': {7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
                  18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29
                  },
           'z': {2.87, 4.08, 3.51, 3.26, 2.51, 2.24, 2.18, 2.07,
                 2.03, 1.67, 1.38, 1.36, 1.19, 1.02, 0.89, 0.78,
                 0.72, 0.64, 0.56, 0.53, 0.48, 0.46, 0.45
                 },
           'dens': {0.738, 0.765, 0.781, 0.792, 0.796, 0.810, 0.825,
                    0.836, 0.842, 0.849, 0.845, 0.848, 0.858, 0.863,
                    0.868, 0.873, 0.877, 0.881, 0.885, 0.889, 0.893,
                    0.897, 0.900
                    },
           'mw': {96, 107, 121, 134, 147, 161, 175, 190, 206, 222,
                  237, 251, 263, 275, 291, 305, 318, 331, 345, 359,
                  374, 388, 402, 449
                  }
           }

    def __init__(self):
        ...

    def intro_fit(self, x: str,
                  y: str):
        """Call function
        Parameters
        ----------
        x, y: string
            default property of the fluid required to adjust parameters
            It can be: cn, z, dens, mw
        """
        set_x = self.com[x]
        set_y = self.com[y]

        x, y = list(set_x), list(set_y)
        return x, y

    def carbon_number_foo(self, z, A=None, B=None):
        """Carbon number fraction function

        Estimates the carbon number based on a composition
        and two fit  parameters
        Use the Pedersen's ratio
        This means that it applies from carbon 6 onward

        Parameters
        ----------
        z: set
          Molar fraction of component
        A, B: float
            A, B fit parameters

        Returns
        -------
        carbon_number: float
            number of single carbon number

        """
        foo_fit = Foo_fit()
        foos = Foos()

        if A is None:
            cn_adj, z_adj = foos.intro_fit('cn', 'z')
            A = foo_fit.fit_AB(cn_adj, z_adj)[0]

        if B is None:
            cn_adj, z_adj = foos.intro_fit('cn', 'z')
            B = foo_fit.fit_AB(cn_adj, z_adj)[1]

        return A + (B*np.log(z))


class Distribution_pedersen(Foos):
    def __init__(self):
        ...

    def p_density(self, cn, L=None, M=None):
        """ Densities function

        Estimates the densities from C6 onwards
        based on increase with carbon number

        Parameters
        ----------
        cn: set
            Carbon number function
        L, D: float
            fit parameters

        Returns
        -------
        density: set
            Density of a given carbon number fraction

        """

        foo_fit = Foo_fit()
        if L is None and M is None:
            L, M = Foos().intro_fit(foo_fit.fit_LM, 'cn', 'dens')

        return L + M * np.log(cn)

    def p_molar_fraction(self, cn, A=None, B=None):
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

        if A is None and B is None:
            A, B = self.A, self.B

        return np.exp((cn-A)/B)

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

        return cn * 14 - 4


class Distribution_cismondi(Foos):
    def __init__(self):
        ...

    def c_molar_fraction(self, cn, Ac=None, Bc=None):
        """ Molar fraction function

        Estimates the molar fraction based on carbon number
        and two fit  parameters

        Parameters
        ----------
        cn: float or set
          Carbon number
        Ac, Bc: float
          fit parameters

        Returns
        -------
        molar_fraction: float
          Molar fraction of a given carbon number fraction
        """
        foo_fit = Foo_fit()

        if Ac is None:
            Ac = Foos().intro_fit(foo_fit.fit_AcBc,
                                  'cn', 'z')[0]
        if Bc is None:
            Bc = Foos().intro_fit(foo_fit.fit_AcBc,
                                  'cn', 'z')[1]

        return np.exp(Ac*cn + Bc)

    def c_molecular_weight(self, cn, C=None):

        """ Molecular weight function

        Estimates the molecular weight based on carbon number
        and the third fit parameter of Cismondi  "C"

        Parameters
        ----------
        cn: float or set
            Carbon number function
        C: float
            fit parameter (third parameter of cismondi et al.)

        Returns
        -------
        molecular_weight: float
            Molecular weight of a given carbon number fraction
        """
        foo_fit = Foo_fit()

        if C is None:
            C = Foos().intro_fit(foo_fit.fit_C, 'cn', 'mw')

        return 84 + C * (cn - 6)

    def c_density(self, cn, Ad=None):
        """ Density function
        Estimates the densities from C6 onwards
        based on [...]

        Parameters
        ----------
        cn: float or set
            Carbon number function
        Ad: float
            fit parameters of cismondi's propose

        Returns
        -------
        density: float
            Density of a given carbon number fraction
        """
        foo_fit = Foo_fit()

        if Ad is None:
            car = self.com['cn']
            d = self.com['dens']
            Ad = foo_fit.fit_Ad(car, d)

        Bd = 0.685 - Ad * np.exp(-0.6)

        return Ad*(np.exp(-cn/10))+Bd


class Foo_fit(Distribution_pedersen, Distribution_cismondi):

    def __init__(self):
        ...

    def fit_AB(self, cn, z):
        """ Parameter setting function

        Parameters
        ----------
        cn: set
            carbon number
        z: set
            molar fraction

        Return
        ------
        A,B: float
            fit parameters
        """

        foos = Foos()
        self.A, self.B = curve_fit(foos.carbon_number_foo,
                                   z, cn,
                                   check_finite=bool)[0]

        return self.A, self.B

    def fit_LM(self, cn, density):
        """ Parameter setting function

        Parameters
        ----------
        cn: set
            carbon number
        density: set
            density

        Return
        ------
        L,M: float
            fit parameters
        """

        distribution_pedersen = Distribution_pedersen()
        self.L, self.M = curve_fit(distribution_pedersen.p_density,
                                   cn, density, check_finite=bool)[0]

        return self.L, self.M

    def fit_C(self, cn, mw):
        """ Parameter setting function
        Parameters
        ----------
        cn: set
            carbon number
        mw: set
            molecular_weight
        Returns
        -------
        C: float
            Third parameter of Cismondi's model
        """

        distribution_cismondi = Distribution_cismondi()
        self.C = curve_fit(distribution_cismondi.c_molecular_weight,
                           cn, mw, check_finite=bool)[0]

        return self.C

    def fit_Ad(self, cn, d):
        """ Parameter setting function
        Parameters
        ----------
        x: set
            carbon number
        d: set
            density for Ad
        Returns
        -------
        C: float
            Third parameter of Cismondi's model
        """
        distribution_cismondi = Distribution_cismondi()
        self.Ad = curve_fit(distribution_cismondi.c_density,
                            cn, d, check_finite=bool)[0]

        return self.Ad

    def fit_AcBc(self, cn, mf):
        """ Parameter setting function
        Parameters
        ----------
        cn: set
            carbon number
        mf: set
            molar_fraction for Ac, Bc

        Returns
        -------
        self.Ac, self.Bc: float
            fit parameter
        """
        distribution_cismondi = Distribution_cismondi()
        self.Ac, self.Bc = curve_fit(distribution_cismondi.c_molar_fraction,
                                     cn, mf, check_finite=bool)[0]

        return self.Ac, self.Bc


class Correlations:
    dic_coeff = {
        "PR": {'c1': 73.4043, 'c2': 97.352,
               'c3': 0.618744, 'c4': -2059.32,
               'd1': 0.0728462, 'd2': 2.18811,
               'd3': 163.91, 'd4': -4043.23, 'd5': 25,
               'e1': 0.373765, 'e2': 0.00549269,
               'e3': 0.0117934, 'e4': -0.00000493049},

        "SRK": {'c1': 163.12, 'c2': 86.052,
                'c3': 0.43375, 'c4': -1877.4,
                'd1': -0.13408, 'd2': 2.5019, 'd3': 208.46,
                'd4': -3987.2, 'd5': 2.0,
                'e1': 0.74310, 'e2': 0.0048122,
                'e3': 0.0096707, 'e4': -0.0000037184}
                  }

    dic_quadratic = {"SRK": {'Ce': 0.480, 'Be': 1.574,
                             'Ae': -0.176},
                     "PR": {'Ce': 0.37464, 'Be': 1.54226,
                            'Ae': -0.26992}}

    def __init__(self, EoS):
        self.dic_coeff = self.dic_coeff[EoS]
        self.dic_quadratic = self.dic_quadratic[EoS]

    def critical_temperature(self, molecular_weight, density,
                             c1=None, c2=None, c3=None, c4=None):
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
        if c1 is None:
            c1 = self.dic_coeff['c1']
        if c2 is None:
            c2 = self.dic_coeff['c2']
        if c3 is None:
            c3 = self.dic_coeff['c3']
        if c4 is None:
            c4 = self.dic_coeff['c4']

        return (
                c1 * density + c2 * (np.log(molecular_weight))
                + c3 * molecular_weight + (c4 / molecular_weight)
                )

    def critical_pression(self, molecular_weight, density,
                          d1=None, d2=None, d3=None,
                          d4=None, d5=None):
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
        if d1 is None:
            d1 = self.dic_coeff['d1']
        if d2 is None:
            d2 = self.dic_coeff['d2']
        if d3 is None:
            d3 = self.dic_coeff['d3']
        if d4 is None:
            d4 = self.dic_coeff['d4']
        if d5 is None:
            d5 = self.dic_coeff['d5']

        pression_log = d1 + d2*(
                density**d5) + (
                d3/molecular_weight) + (
                d4/(molecular_weight**2)
                )

        return np.exp(pression_log)


class Proper_plus():
    def __init__(self):
        ...

    def properties_plus(self, molarfraction_values, molecularweight_values):
        """Function to calculate residual fraction properties

        Parameters
        ----------
        molarfraction_values: float
            Molar fraction values to i
        molecularweight_values: float
            Molecular weight values to i

        Return
        ------
        molecularplus: float
            Molecular values to residual fraction
        molarplus: float
            Molar fraction to residual fraction
        """
        molecularplus = np.array([])
        sa = np.array([])
        mfv = np.array([])

        for i, j in zip(molecularweight_values,
                        molarfraction_values):
            sa = np.append(sa, i*j)
            mfv = np.append(mfv, j)
            molecularplus = np.append(molecularplus,
                                      (sa.sum())/mfv.sum())

        molarfractionplus = np.array([])
        su = np.array([])

        for i in molarfraction_values:
            su = np.append(su, i)
            molarfractionplus = np.append(molarfractionplus,
                                          su.sum())

        return molecularplus, molarfractionplus


class Residual_fraction(Proper_plus, Foo_fit):

    def __init__(self, mw_max, mf_max):
        self.molecularweight_max = mw_max
        self.molarfraction_max = mf_max

    def carbon_number_max(self, carbon_range, A, B, C):
        """Maximum carbon number based on Cismondi's observations

        Parameters
        ----------
        carbon_range: set
            Carbon range to distribute
        Ac: float
            fit parameter
        Bc: float
            fit parameter
        C: float
            fit parameter
        Returns
        -------
        carbonnumber_max: int
            Maximun carbon number
        """

        distribution_cismondi = Distribution_cismondi()
        distribution_pedersen = Distribution_pedersen()

        proper_plus = Proper_plus()

        carbonnumber_max = carbon_range[0]
        molarfraction_values = np.array(
            distribution_pedersen.p_molar_fraction(
                carbon_range, A, B))
        molecularweight_values = np.array(
            distribution_cismondi.c_molecular_weight(
                carbon_range, C))

        molecularplus, molarfractionplus = proper_plus.properties_plus(
            molarfraction_values, molecularweight_values)

        for i, j in zip(molecularplus, molarfractionplus):
            carbonnumber_max += 1

            if i > self.molecularweight_max or j > self.molarfraction_max:
                break

        return carbonnumber_max


class Lumping():
    def __init__(self):
        ...

    def criteriotanto(self, dato):
        ...

    def lumy(self, limits, datos, colum_name):
        """Function to inset ID index to lumpy func.

        Parameters
        ----------
        limits: list of tuples
            ranges limits
        datos: set
            distribution data set
        colum_name: set
            data frame column to realize lumping

        Return
        ------
        lumy_set: set
            data set with ID refers to lumping
        """
        count = -1
        co = np.array([])

        for i in datos[colum_name]:
            count += 1
            if i in range(limits[count][0],
                          limits[count][1]+1,
                          1):
                co = np.append(co, count)
                count -= 1
            else:
                co = np.append(co, count+1)

        datos['ID'] = co

        datos.set_index([pd.Index(co), 'ID'])

        return datos

    def lumpy(self, limits, datos, colum_name):
        """ Function to det. carbon ranges for lumping

        Parameters
        ----------
        limits: tupla
            ranges limits
        datos: set
            distribution data set
        colum_name: set
            data frame column to realize lumping

        Return
        ------
        lumping_set: set
            lumping data set (groupy)

        """
        lumping = Lumping()
        df = lumping.lumy(limits, datos, colum_name)

        data_lumping = df.groupby(by=["ID"],
                                  dropna=False).mean()

        data_lumping[colum_name] = limits

        return data_lumping


def test_introfit_001():
    x = [7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
         18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
         28, 29]

    assert set(Foos().intro_fit('cn', 'z')[0]) == set(x)


def test_introfit_002():
    y = [2.87, 4.08, 3.51, 3.26, 2.51, 2.24, 2.18, 2.07,
         2.03, 1.67, 1.38, 1.36, 1.19, 1.02, 0.89, 0.78,
         0.72, 0.64, 0.56, 0.53, 0.48, 0.46, 0.45]

    assert set(Foos().intro_fit('cn', 'z')[1]) == set(y)


def test_cn():
    z = np.array([2.87])
    assert Foos().carbon_number_foo(z) == 7