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


def f(x, A, B):
    return np.exp(-A*x)+np.exp(-B*x)



# Messwerte #########################################################

t_V, N_V = np.genfromtxt('Vanadium.dat' , unpack=True)
t_Rh, N_Rh = np.genfromtxt('Rhodium.dat', unpack=True)

N_U_V= np.array([129,129,129,129,129,129,129,129,129,129,143,143,143,143,143,143,143,143,143,143,144,144,144,144,144,144,144,144,144,144,136,136,136,136,136,136,136,136,136,136,139,139,139,139])
N_U_V_err = unp.uarray(N_U_V, np.sqrt(N_U_V))

N_V_err = unp.uarray(N_V, np.sqrt(N_V))

N_V_ohne_U = N_V_err-N_U_V_err/10



N_U_Rh =np.array([129,129,129,129,129,129,129,129,129,129,129,129,129,129,129,129,129,129,129,129,143,143,143,143,143,143,143,143,143,143,143,143,143,143,143,143,143,143,143,143,136,136,136,136])
N_U_Rh_err =unp.uarray(N_U_Rh,np.sqrt(N_U_Rh))

N_Rh_err = unp.uarray(N_Rh, np.sqrt(N_Rh))

N_Rh_ohne_U = N_Rh_err-N_U_Rh_err/20




# Plots #####################################################################

# Plot von Vanadium %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par, cov= np.polyfit(t_V, np.log(noms(N_V_ohne_U)), deg=1, cov=True)
err= np.sqrt(np.diag(cov))

plt.errorbar(t_V, noms(N_V_ohne_U),xerr=stds(N_V_ohne_U) ,fmt='ko', label='Messwerte')
x_plot = np.linspace(0, 1400, 10000)
plt.plot(x_plot, np.exp(x_plot*par[0]+par[1]), 'r-', label='Fitkurve')
plt.yscale('log')
plt.legend(loc="best")
plt.xlabel(r' $t \:/\:s$')
plt.ylabel(r' $I\:/\:Imp/30s$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_V.pdf')
plt.close()

par_V = unp.uarray(par, err)

# Plot von Rhodium %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print(noms(t_Rh[26:]))

par, cov = np.polyfit( t_Rh[22:], np.log(noms(N_Rh_ohne_U[22:])), deg=1, cov=True)

par2, cov2 = curve_fit(
    f,
    t_Rh,
    np.log(noms(N_Rh_ohne_U)),
    sigma=None,
    absolute_sigma=True,
    p0=[100, 1e-02]
)



plt.errorbar(t_Rh, noms(N_Rh_ohne_U),xerr=stds(N_Rh_ohne_U) ,fmt='kx', label='Messwerte')
x_plot = np.linspace(0, 660, 10000)
plt.plot(x_plot, np.exp(x_plot*par[0]+par[1]) , 'r-', label='Gerade des langlebigen Zerfalls')
#plt.plot(x_plot, np.exp(-par2[0]*x_plot)+np.exp(-par2[1]*x_plot), 'b-', label='Fitkurve' )
plt.yscale('log')
plt.legend(loc="best")
plt.xlabel(r' $t \:/\:s$')
plt.ylabel(r' $I\:/\:Imp/15s$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_Rh.pdf')
plt.close()
