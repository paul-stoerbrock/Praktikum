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

def c(E0, E, m0, c, tau):
    return np.sqrt(E0**2-E**2)/(m0* c) *tau

def d(a, b, E, E0):
    return - 1/b * np.log((a+b*E)/(a+b*E0))


#with open('build/luft.tex', 'w') as f:
#  f.write(make_SI(l*1e-06,r'\giga\electronvolt\per\meter', figures=2))


x_plot = np.linspace(0, 140000000, 10000)
plt.plot(x_plot, d( 2.302, 3.617*1e-6, 105, x_plot), 'r-', label='Energieverlust')
plt.legend(loc="best")
plt.xlabel(r' Energie in MeV')
plt.ylabel(r' Reichweite in cm')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_d.pdf')
plt.close()


x_plot = np.linspace(105, 400, 10000)
plt.plot(x_plot, c(x_plot, 105, 1.8*1e-28, noms(const.c*100), 2.19*1e-6 ), 'r-', label='Zerfall')
plt.legend(loc="best")
plt.xlabel(r' Energie in MeV')
plt.ylabel(r' Reichweite in cm')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_c.pdf')
plt.close()