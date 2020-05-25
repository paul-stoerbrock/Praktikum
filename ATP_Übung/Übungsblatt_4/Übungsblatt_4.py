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
from astropy import constants as astro
import scipy.constants as const #Bsp.: const.physical_constants["proton mass"], output -> value, unit, error
import math


def make_SI(num, unit, exp='', figures=None):
    ''' Format an uncertainties ufloat as a \SI quantity '''
    if np.any(stds([num])):
        if figures is None:
            figures = ''
        x = '{0:.{1:}uf}'.format(num, figures).replace('/', '')
    else:
        x = '{0:.{1:}f}'.format(num, figures)

    return r'\SI{{{}{}}}{{{}}}'.format(x, exp, unit)

M = -26.74-5*math.log10(149.6*10**6/(10*3.086*1e+13))

r = 10**((-1.46-1.43)/5+1)
print(r)

# tex file for M

with open('build/M.tex', 'w') as f:
  f.write(make_SI(M,r'', figures=4))

# tex file for r

with open('build/r.tex', 'w') as f:
  f.write(make_SI(r,r'pc', figures=4))