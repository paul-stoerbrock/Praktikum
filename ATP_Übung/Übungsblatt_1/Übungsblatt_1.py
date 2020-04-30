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

#def A_S(R_S):                       #Sonnenoberfläche
#    return 4*np.pi*R_S**2
#
#def A_E(R_E):                       #Erdoberfläche
#    return 4*np.pi*R_E**2
#
#def P_S(A_S, T_S):                  #Luminosität (Sonne)
#    return const.Stefan_Boltzmann*A_S*T_S**4

def L(P_S, R):                      #Flussdichte im Abstand 1 au
    return P_S/(4*np.pi*R)

T_S = 5780
R_S = 6.69*10**5
R_E = 6360
R   = const.au

A_S = 4*np.pi*R_S**2
A_E = 4*np.pi*R_E**2
P_S = const.sigma*A_S*T_S**4

# tex file of T_S

with open('build/T_S.tex', 'w') as f: 
  f.write(make_SI(T_S, r'\kelvin' ,figures=0))

# tex file of R_S

with open('build/R_S.tex', 'w') as f: 
  f.write(make_SI(R_S, r'\kilo\meter' ,figures=0))

# tex file of R_E

with open('build/R_E.tex', 'w') as f: 
  f.write(make_SI(R_E, r'\kilo\meter' ,figures=0))

# tex file of A_S

with open('build/A_S.tex', 'w') as f: 
  f.write(make_SI(A_S, r'\kilo\meter\tothe{2}' ,figures=0))

# tex file of A_E

with open('build/A_E.tex', 'w') as f: 
  f.write(make_SI(A_E, r'\kilo\meter\tothe{2}' ,figures=0))

# tex file of P_S

with open('build/P_S.tex', 'w') as f: 
  f.write(make_SI(P_S*1e-15, r'\peta\watt' ,figures=1))