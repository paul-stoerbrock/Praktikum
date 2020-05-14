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

def make_SI(num, unit, exp='', figures=None):
    ''' Format an uncertainties ufloat as a \SI quantity '''
    if np.any(stds([num])):
        if figures is None:
            figures = ''
        x = '{0:.{1:}uf}'.format(num, figures).replace('/', '')
    else:
        x = '{0:.{1:}f}'.format(num, figures)

    return r'\SI{{{}{}}}{{{}}}'.format(x, exp, unit)





# = np.genfromtxt('Bragg.dat', unpack = True)
# = np.genfromtxt('Brom.dat', unpack = True)
# = np.genfromtxt('Emissionsspektrum.dat', unpack = True)
# = np.genfromtxt('Gallium.dat', unpack = True)
# = np.genfromtxt('Rubidium.dat', unpack = True)
# = np.genfromtxt('Strontium.dat', unpack = True)
# = np.genfromtxt('Zink.dat', unpack = True)
# = np.genfromtxt('Zirkonium.dat', unpack = True)
