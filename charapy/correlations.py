import numpy as np

coefficients = {
                'PR':{'c1': 73.4043, 'c2': 97.352,
                      'c3': 0.618744, 'c4': -205932,
                      'd1': 0.0728462, 'd2': 2.18811,
                      'd3': 163.91, 'd4': -4043.23, 'd5': 25,
                      'e1': 0.373765, 'e2': 0.00549269,
                      'e3': 0.0117934, 'e4': -0.00000493049},
                'SRK':{'c1': 163.12, 'c2': 86.052,
                       'c3': 0.43375, 'c4': -1877.4,
                       'd1': -0.13408, 'd2': 2.5019, 'd3':208.46,
                       'e1': 0.74310, 'e2': 0.0048122,
                       'e3': 0.0096707, 'e4': -0.0000037184}
                }

quadratic_coef = {
                  'SRK': {'Ce': 0.480, 'Be': 1.574,
                          'Ae': -0.176},
                  'PR': {'Ce': 0.37464, 'Be': 1.54226,
                         'Ae': -0.26992}
                  }

def __init__(self,EoS):
    self.dic_coeff = self.coefficients[EoS]
    self.quadratic_coef = self.quadratic_coef[EoS]

def critical_temperature(self, mw, rho,
                         c1=None, c2=None,
                         c3=None, c4=None):

    '''Critical temperature correlation

    Parameters
    ----------
    mw: float
        molecular weight
    rho: float
        density
    c1, c2, c3, c4: float
        parameters

    Returns
    -------
    cT = float
        critical temperature
    '''

    if c1 is None: c1 = coefficients['c1']
    if c2 is None: c2 = coefficients['c2']
    if c3 is None: c3 = coefficients['c3']
    if c4 is None: c4 = coefficients['c4']

    cT = (c1*rho+c2*(np.log(mw))
          +c3*mw+(c4/mw))

    return cT

def critical_presion(self,mw,rho,
                     d1=None, d2=None, d3=None,
                     d4=None, d5=None):
    '''Critical pression correlation

    Parameters
    ----------
    mw: float
        molecular weight
    rho: float
        density
    d1, d2, d3, d4, d5: float
        fit parameters

    Returns
    -------
    cP = float
        critical pression
    '''

    if d1 is None: d1 = coefficients['d1']
    if d2 is None: d2 = coefficients['d2']
    if d3 is None: d3 = coefficients['d3']
    if d4 is None: d4 = coefficients['d4']
    if d5 is None: d5 = coefficients['d5']

    p_log = d1+d2*(rho**d45)+(
            d3/mw)+(d4/(mw**2))

    cP = np.exp(p_log)
    return cP