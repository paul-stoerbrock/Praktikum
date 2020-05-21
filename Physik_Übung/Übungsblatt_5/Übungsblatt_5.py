import pylab
import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.unumpy import (
    nominal_values as noms,
    std_devs as stds,
)
from scipy import stats 
from scipy.stats import sem #standard error of mean = sem(x)
from scipy.optimize import curve_fit #function curve_fit 
import scipy.constants as const #Bsp.: const.physical_constants["proton mass"], output -> value, unit, error

# Funktionsdefinitionen ########################################################################################################################################################

def make_SI(num, unit, exp='', figures=None):
    ''' Format an uncertainties ufloat as a \SI quantity '''
    if np.any(stds([num])):
        if figures is None:
            figures = ''
        x = '{0:.{1:}uf}'.format(num, figures).replace('/', '')
    else:
        x = '{0:.{1:}f}'.format(num , figures)

    return r'\SI{{{}{}}}{{{}}}'.format(x, exp, unit)


P_1 = np.exp(-4 * np.sqrt(-2 *const.m_e *(-4.5*const.e)**3)/(3*const.hbar*5e08*const.e))
P_2 = np.exp(-4 * np.sqrt(-2 *const.m_e *(-4.5*const.e)**3)/(3*const.hbar*5e09*const.e))

# tex file for P_1

with open('build/P_1.tex', 'w') as f:
  f.write(make_SI(P_1*1e+57,r'\percent', exp='1e-55', figures=1))

# tex file for P_2

with open('build/P_2.tex', 'w') as f:
  f.write(make_SI(P_2*1e+06,r'\percent', exp='1e-04', figures=1))



print(P_1)
print(P_2)
