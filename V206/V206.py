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
        x = '{0:.{1:}f}'.format(num, figures)

    return r'\SI{{{}{}}}{{{}}}'.format(x, exp, unit)




# Definition der Messdaten ##############################################################################################

t, pa, T2, pb, T1, P = np.genfromtxt('data.txt', unpack=True)

t_s =t*60
pa_Pa =pa*1e05+1e05
pb_Pa = pb*1e05+1e05
T2_K = const.convert_temperature(T2, 'c', 'K')
T1_K = const.convert_temperature(T1, 'c', 'K')


# Erstellung der Plots ######################################################################################################

# Plot für T1(t) #############################################################################################
par, covm = np.polyfit(t_s, T1_K, deg=2, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(t_s, T1_K,'kx', label='Messwerte')
x_plot = np.linspace(0, 1800, 5000)
plt.plot(x_plot, par[0]*x_plot**2 + par[1]*x_plot + par[2], 'r-', label="Lineare Regression")
plt.legend(loc="best")
plt.xlabel(r'Zeit $t/s$')
plt.ylabel(r'Temperatur $T/K$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotT1.pdf')
plt.close()

parT1=unp.uarray(par, err)

print(T2_K)
# Plot für T2(t) #############################################################################################
par, covm = np.polyfit(t_s, T2_K, deg=2, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(t_s, T2_K,'kx', label='Messwerte')
x_plot = np.linspace(0, 1800, 5000)
plt.plot(x_plot, par[0]*x_plot**2+par[1]*x_plot+par[2], 'r-', label="Lineare Regression")
plt.legend(loc="best")
plt.xlabel(r'Zeit $t/s$')
plt.ylabel(r'Temperatur $T/K$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotT2.pdf')
plt.close()

parT2=unp.uarray(par, err)

