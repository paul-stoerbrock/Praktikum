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
import math
import sympy

def make_SI(num, unit, exp='', figures=None):
    ''' Format an uncertainties ufloat as a \SI quantity '''
    if np.any(stds([num])):
        if figures is None:
            figures = ''
        x = '{0:.{1:}uf}'.format(num, figures).replace('/', '')
    else:
        x = '{0:.{1:}f}'.format(num, figures)

    return r'\SI{{{}{}}}{{{}}}'.format(x, exp, unit)


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
        
    return latex(sympy.sqrt(s), symbol_names=latex_names) # automatic Gauß error function

def rel_err(mess, lit):
    return abs((unp.nominal_values(mess)-unp.nominal_values(lit))/unp.nominal_values(lit))*100



# Messwerte ######################################################################

UA_T1, IA_UA_T1 = np.genfromtxt('dataUA_T1.txt' ,unpack = True) # UA_T1 bei 25,3°C gemessen bei UB =11V, UA_T2 bei 145,5°C
UA_T2, IA_UA_T2 = np.genfromtxt('dataUA_T2.txt' ,unpack = True) # UA_T2 bei 145,5°C
UB_T1 , IA_UB_T1, UB_T2, IA_UB_T2 = np.genfromtxt('dataUB.txt' ,unpack = True) # UB_T1 bei 164°C und UA=1V , UB_T2 bei 175°C
 
T1_UA = const.convert_temperature( 25.3, 'C' ,'K')
T2_UA = const.convert_temperature( 145.5, 'C' ,'K')
U_B = 11

T1_UB = const.convert_temperature( 164, 'C' ,'K')
T2_UB = const.convert_temperature( 175, 'C' ,'K')
U_A = 1




# Plots ####################################################################


# Plot für die Steigung der Gegenspannung %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par1, cov1 = np.polyfit( UA_T1, IA_UA_T1, deg=1, cov=True )

err1 = np.sqrt(np.diag(cov1))

par2, cov2 = np.polyfit( UA_T2, IA_UA_T2, deg=1, cov=True )

err2 = np.sqrt(np.diag(cov2))


plt.plot(UA_T1, IA_UA_T1 ,'kx' , label='Messwerte bei 25,3°C')
plt.plot(UA_T2, IA_UA_T2 ,'bx' , label='Messwerte bei 145,5°C')
x_plot = np.linspace(0, 5, 1000)
plt.plot(x_plot, par1[0]*x_plot+par1[1], 'r-', label='Fitgerade bei 25,3°C')
plt.plot(x_plot, par2[0]*x_plot+par2[1], 'b-', label='Fitgerade bei 145,5°C')
plt.legend(loc="best")
plt.xlabel(r' $U_A \:/\:V$')
plt.ylabel(r'$I_A \:/\: nA$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_UA.pdf')
plt.close()


par_T1 = unp.uarray(par1, err1)

par_T2 = unp.uarray(par2, err2)



# Plot für Die Beschleunigungsspannung %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plt.plot(UB_T1, IA_UB_T1 ,'kx' , label='Messwerte bei 164°C')
plt.plot(UB_T2, IA_UB_T2 ,'bx' , label='Messwerte bei 175°C')
plt.legend(loc="best")
plt.xlabel(r' $U_B \:/\:V$')
plt.ylabel(r'$I_A \:/\: nA$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_UB.pdf')
plt.close()
