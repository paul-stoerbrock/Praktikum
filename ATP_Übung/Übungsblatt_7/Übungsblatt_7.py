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


def f(n):
    return (1/(np.sqrt(1-1/n**2))) *const.c**2 *const.m_e



l=f(1.0003)/const.e
w=f(1.333)/const.e
# tex file for Luft

with open('build/luft.tex', 'w') as f:
  f.write(make_SI(l*1e-06,r'\giga\electronvolt\per\meter', figures=2))

with open('build/wasser.tex', 'w') as f:
  f.write(make_SI(w*1e-06,r'\giga\electronvolt\per\meter', figures=2))