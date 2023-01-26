import numpy as np
from charapy.distributions import cismondi, pedersen

class Proper_plus():
    def __init__(self):
        ...

    def prop_plus(self,z_values,mw_values):
        '''Function to calculate residual
        fraction properties

        Parameters
        ----------
        z_values: float
            molar fraction values
        mw_values: float
            molecular weight values

        Return
        ------
        molecularplus: float
            molecular weight values
            to residual fraction
        molarplus: float
            molar fraction to residual fraction
        '''

        mwplus = np.array([])
        sa = np.array([])
        mfv = np.array([])
        zplus = np.array([])
        su = np.array([])

        for i, j in zip(mw_values, z_values):
            sa = np.append(sa, i*j)
            mfv = np.append(mfv, j)
            mwplus = np.append(mwplus,sa.sum()/mfv.sum())

        for i in z_values:
            su = np.append(su, i)
            zplus = np.append(zplus,su.sum())

        return mwplus, zplus

class Residual_fraction(Proper_plus):

    def __init__(self, mw_max, z_max):
        self.mw_max = mw_max
        self.z_max = z_max

    def cn_max(self, cn_range, A, B,fun,*args):
        '''Function to calculated max carbon number
        based on Cismondi's observations

        Parameters
        ----------
        cn_range: int
            carbon number range to distribute
        A, B, C: float
            fit parameters

        Returns
        -------
        cn_max: int
            max carbon number
        '''

        dcis = cismondi.Distribution_cismondi()
        dped = pedersen.Distribution_pedersen()

        cn_max = cn_range[0]
        C = args

        if fun == 'cis':
            mw_values = np.array(dcis.molecular_weight(
                             cn_range,C))
            z_values = np.array(dcis.molar_fraction(
                            cn_range,A,B))

        elif fun == 'ped':
            mw_values = np.array(dped.molecular_weight(cn_range))
            z_values = np.array(dped.molar_fraction(
                            cn_range,A,B))

        else: raise SyntaxError('fun puede ser: ped o cis')

        pplus = Proper_plus()
        mwplus, zplus = pplus.prop_plus(z_values,mw_values)

        for i, j in zip(mwplus, zplus):
            cn_max += 1

            if round(i,5) >= self.mw_max or round(j,5) >= self.z_max:
                break

        return cn_max