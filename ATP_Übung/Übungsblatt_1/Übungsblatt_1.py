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
R_S = 6.96*10**8
R_E = 6360*1e03
R   = const.au

A_S = 4*np.pi*R_S**2
A_E = 4*np.pi*R_E**2
P_S = const.sigma*A_S*T_S**4
Lum = L(P_S,R) 


# Zu Aufgabe 3

f_H1 = 1/(21.1*1e-02) * const.c
E_H1 = f_H1*2*np.pi*const.physical_constants["Planck constant over 2 pi"][0]

print(P_S)
print(E_H1 /const.e )

# tex file of T_S

with open('build/T_S.tex', 'w') as f: 
  f.write(make_SI(T_S, r'\kelvin' ,figures=0))

# tex file of R_S

with open('build/R_S.tex', 'w') as f: 
  f.write(make_SI(R_S*1e-03, r'\kilo\meter' ,figures=0))

# tex file of R_E

with open('build/R_E.tex', 'w') as f: 
  f.write(make_SI(R_E*1e-03, r'\kilo\meter' ,figures=0))

# tex file of A_S

with open('build/A_S.tex', 'w') as f: 
  f.write(make_SI(A_S, r'\kilo\meter\tothe{2}' ,figures=0))

# tex file of A_E

with open('build/A_E.tex', 'w') as f: 
  f.write(make_SI(A_E, r'\kilo\meter\tothe{2}' ,figures=0))

# tex file of P_S

with open('build/P_S.tex', 'w') as f: 
  f.write(make_SI(P_S*1e-26, r'\watt',exp='e26' ,figures=2))

# tex file of F_H1

with open('build/f_H1.tex', 'w') as f: 
  f.write(make_SI(f_H1*1e-09,r'\giga\hertz' ,figures=2))


# tex file of E_H1

with open('build/E_H1.tex', 'w') as f: 
  f.write(make_SI(E_H1, r'\electronvolt' ,figures=2))
