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


def rel_err(mess, lit):
    return abs((unp.nominal_values(mess)-unp.nominal_values(lit))/unp.nominal_values(lit))*100



# Messdaten #################################################################

U_kenn, N = np.genfromtxt('Kennlinie.dat', unpack=True)

U_zaehl, I = np.genfromtxt('Zaehlrohrstrom.dat', unpack=True)


N_err = np.sqrt(N)

N_miterr = unp.uarray(N, N_err)

I_miterr = unp.uarray(I*1e-06, 0.05e-06)

# Plots ######################################################################

# Plot zur Kennlinie des Geiger-Müller Zählrohrs

plt.errorbar(U_kenn, N, yerr=stds(N_miterr) ,fmt='kx', label='Messwerte mit Fehler')
x_plot = np.linspace(26, 30, 1000)
#plt.plot(x_plot, , 'r-', label='Fitkurve')
plt.legend(loc="best")
plt.xlabel(r'Spannung $U \:/\:V$')
plt.ylabel(r'Intensität $I\:/\:Imp/s$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_kenn.pdf')
plt.close()

