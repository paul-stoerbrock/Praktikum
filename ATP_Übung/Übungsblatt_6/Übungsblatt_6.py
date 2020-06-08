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

h, I = np.genfromtxt('data.txt', unpack=True)








def I_h(x, a, b, c):
  return a*x+b*x**2+c*x**3

popT, pcov = curve_fit(
    I_h,
    h,
    I,
    sigma=None,
    absolute_sigma=True,
    p0=[100, 28, -30, 20]
    )
err_var = np.sqrt(np.diag(pcov))


plt.plot(h, I, 'kx', label='Messwerte')
x_plot = np.linspace(0, 100, 1000)
plt.plot(x_plot, Gauß(x_plot,*popT), linestyle='-', label='Ausgleichskurve')

plt.xlabel(r'I')
plt.ylabel(r'h')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plot18.pdf')
plt.close()

