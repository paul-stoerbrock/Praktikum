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
import sympy

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

def error(f, err_vars=None):
    from sympy import Symbol, latex
    s = 0
    latex_names = dict()
    
    if err_vars == None:
        err_vars = f.free_symbols
        
    for v in err_vars:
        err = Symbol('latex_std_' + v.name)
        s += f.diff(v)**2 * err**2
        latex_names[err] = '\\sigma_{' + latex(v) + '}'
        
    return latex(sympy.sqrt(s), symbol_names=latex_names)



# Messwerte #########################################################

t_V, N_V = np.genfromtxt('Vanadium.dat' , unpack=True)
t_Rh, N_Rh = np.genfromtxt('Rhodium.dat', unpack=True)


N_V_err = unp.uarray(N_V, np.sqrt(N_V))
N_Rh_err = unp.uarray(N_Rh, np.sqrt(N_Rh))




# Plots #####################################################################

# Plot von Vanadium %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plt.errorbar(t_V, N_V,xerr=stds(N_V_err) ,fmt='ko', label='Messwerte')
#x_plot = np.linspace(26, 30, 1000)
#plt.plot(x_plot, gauß(x_plot, *popBragg), 'r-', label='Fitkurve')
plt.legend(loc="best")
plt.xlabel(r'Winkel $t \:/\:s$')
plt.ylabel(r' $I\:/\:Imp/s$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_V.pdf')
plt.close()


# Plot von Rhodium %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plt.errorbar(t_Rh, N_Rh,xerr=stds(N_Rh_err) ,fmt='ko', label='Messwerte')
#x_plot = np.linspace(26, 30, 1000)
#plt.plot(x_plot, gauß(x_plot, *popBragg), 'r-', label='Fitkurve')
plt.legend(loc="best")
plt.xlabel(r'Winkel $t \:/\:s$')
plt.ylabel(r' $I\:/\:Imp/s$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_Rh.pdf')
plt.close()