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


def make_SI(num, unit, exp='', figures=None):
    ''' Format an uncertainties ufloat as a \SI quantity '''
    if np.any(stds([num])):
        if figures is None:
            figures = ''
        x = '{0:.{1:}uf}'.format(num, figures).replace('/', '')
    else:
        x = '{0:.{1:}f}'.format(num, figures)

    return r'\SI{{{}{}}}{{{}}}'.format(x, exp, unit)

R_S = 6.96*10**8                    #Radius in m
A_S = 4*np.pi*R_S**2
T_S = 5780                          #Temperatur in K
P_S = const.sigma*A_S*T_S**4

M_stern = 0.083*1.98840987*1e+30

#Eddington Leuchtkraft:
L = const.G*4*np.pi*M_stern*const.c/0.02

P_sol = P_S/L

# tex file for L

with open('build/L.tex', 'w') as f:
  f.write(make_SI(L*1e-30,r'\watt', figures=2))

# tex file for P_sol

with open('build/P_sol.tex', 'w') as f:
  f.write(make_SI(P_sol,r'', figures=5))

print(astro.M_sun)