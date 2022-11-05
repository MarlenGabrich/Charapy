import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

class Distribution_pedersen():
    """El radio de estimación de Pedersen aplica
    desde Carbono 6 en adelante 
    """

    def __init__(self):
        ...

    def carbon_number_foo(self, z, A, B):
        """Función general para estimar el número de carbono

        Parameters
        ----------
        z: set or float
          fracción molar de los componentes
        A, B: float
            A, B fit parameters

        Returns
        -------
        carbon_number: float or set
            número de carbono

        """

        carbon_number = B*np.log10(z) + A

        return carbon_number.astype(int)

    def p_density(self, cn, L, M):
        """ Función para calcular la densidad

        Parameters
        ----------
        cn: float or set
            número de carbono
        L, D: float
            parámetros de ajuste

        Returns
        -------
        density: float or set
            Densidad del número de carbono dado
        """

        density = L*np.log10(cn) + M

        return density.astype(float)

    def p_molar_fraction(self, cn, A, B):
        """Función para calcular la fracción molar

        Parameters
        ----------
        cn: float or set
            número de carbono
        A,B: float
            parámetros de ajuste

        Returns
        -------
        molar_fraction: float or set
            fracción molar para cada componente
        """

        return np.exp((cn-A)/B)

    def p_molecular_weight(self, cn):
        """ Función para calcular el peso molecular 
        del componente

        Parameters
        ----------
        cn: float or set
            número de carbono

        Returns
        -------
        molecular_weight: float or set
            peso molecular

        """

        return cn * 14 - 4

class Distribution_cismondi():
    def __init__(self):
        ...

    def c_molar_fraction(self, cn, Ac, Bc):
        """ Función para calcular la fracción molar

        Parameters
        ----------
        cn: float or set
          número de carbono
        Ac, Bc: float
          parámetros de ajuste

        Returns
        -------
        molar_fraction: set or float
          fracción molar
        """

        return np.exp(Ac*cn + Bc)

    def c_molecular_weight(self, cn, C):
        """ Función para calcular el peso molecular

        Parameters
        ----------
        cn: float or set
            número de carbono
        C: float
            tercer parámetro de ajuste de Cismondi et al.

        Returns
        -------
        molecular_weight: set or float
            peso molecular
        """
        
        return 84 + C*(cn-6)

    def c_density(self, cn, Ad):
        """ Función para calcular la densidad

        Parameters
        ----------
        cn: float or set
            número de carbono
        Ad: float
            parámetro de ajuste propuesto por Cismondi et al

        Returns
        -------
        density: float or set
            densidad
        """

        Bd = 0.685 - Ad * np.exp(-0.6)

        return Ad*(np.exp(-cn/10))+Bd


class Foo_fit(Distribution_pedersen, Distribution_cismondi):

    def __init__(self):
        ...

    def fit_AB(self, cn_fit, z_fit):
        """ Función de ajuste A, B

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
        x_data = np.log10(z_fit)
        y_data = cn_fit

        A,B = np.polyfit(x_data, y_data,1)

        return A,B

    def fit_LM(self, cn_fit, density_fit):
        """ Función de ajuste L, M

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
        x_data = np.log10(cn_fit)
        y_data = density_fit
        L, M = np.polyfit(x_data, y_data, 1)

        return L, M

    def fit_C(self, cn, mw):
        """ Función de ajuste C

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

        x_data = cn
        y_data = mw

        distribution_cismondi = Distribution_cismondi()
        C = curve_fit(distribution_cismondi.c_molecular_weight,
                           x_data, y_data)[0]

        return C

    def fit_Ad(self, cn, d):
        """ FUnción de ajuste Ad

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
        x_data = cn
        y_data = d

        distribution_cismondi = Distribution_cismondi()
        Ad = curve_fit(distribution_cismondi.c_density,
                            x_data, y_data)[0]

        return Ad

    def fit_AcBc(self, cn, mf):
        """ Función de ajuste Ac, Bc
        
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
        x_data = cn
        y_data = mf

        distribution_cismondi = Distribution_cismondi()
        Ac, Bc = curve_fit(distribution_cismondi.c_molar_fraction,
                                     x_data, y_data)[0]

        return Ac, Bc


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
    y = [2.87, 4.08, 3.51, 3.26, 2.51, 2.24, 2.18, 2.07,
         2.03, 1.67, 1.38, 1.36, 1.19, 1.02, 0.89, 0.78,
         0.72, 0.64, 0.56, 0.53, 0.48, 0.46, 0.45]

    assert 4

def test_introfit_002():
    y = [2.87, 4.08, 3.51, 3.26, 2.51, 2.24, 2.18, 2.07,
         2.03, 1.67, 1.38, 1.36, 1.19, 1.02, 0.89, 0.78,
         0.72, 0.64, 0.56, 0.53, 0.48, 0.46, 0.45]

    assert 4#set(Foos().intro_fit('cn', 'z')[1]) == set(y)

def test_cn():
    assert 4#Foos().carbon_number_foo(2.87) == 7

def test_density_pedersen():
    assert 0.6 <= Distribution_pedersen().p_density(7) <= 0.876

def test_density_cismondi():
    assert 0.6 <= Distribution_cismondi().c_density(7) <= 0.876

def test_molarfraction_pedersen():
    assert 0.95 <= Distribution_pedersen().p_molar_fraction(7) <= 4.79

def test_molarfraction_cismondi():
    assert 0.96 <= Distribution_cismondi().c_molar_fraction(7) <= 4.79

def test_molecular_weight_pedersen():
    assert 94 <= Distribution_pedersen().p_molecular_weight(7) <= 98

def test_molecular_weight_cismondi(): 
    assert 94 <= Distribution_cismondi().c_molecular_weight(7) <= 98