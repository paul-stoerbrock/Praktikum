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

#13a)--------------------------------------------------------------------------------------------------------------------------

Lohbusch_B =342 
Lohbusch_L =114 
Lohbusch_A = Lohbusch_B*Lohbusch_L
N = Lohbusch_A/25
sigma = 0.25
Baumdichte = N/Lohbusch_A

L = 1/(Baumdichte*0.25)

# tex file Länge

with open('build/Länge.tex', 'w') as f:
  f.write(make_SI(Lohbusch_L,r'\meter', figures=0))

# tex file Breite

with open('build/Breite.tex', 'w') as f:
  f.write(make_SI(Lohbusch_B,r'\meter', figures=0))

# tex file Fläche

with open('build/Fläche.tex', 'w') as f:
  f.write(make_SI(Lohbusch_A,r'\meter\squared', figures=0))

# tex file Baumzahl

with open('build/Baumzahl.tex', 'w') as f:
  f.write(make_SI(N,r'\meter', figures=1))

# tex file Baumdichte

with open('build/Baumdichte.tex', 'w') as f:
  f.write(make_SI(Baumdichte,r'\per\meter\squared', figures=2))

# tex file lambda

with open('build/lambda.tex', 'w') as f:
  f.write(make_SI(L,r'\meter', figures=0))

print(Baumdichte)
print(Lohbusch_A)
print(N)
print(L)

#13b)--------------------------------------------------------------------------------------------------------------------------

rho = 10**(-10) *(3.0856*10**16)**(-3)
sigma_13b = np.pi*(6.96342*10**8)**2

L_13b = 1/(rho*sigma_13b)
L_13b_parsec = L_13b*3.2407*10**(-17)

print(L_13b)
print(L_13b_parsec)

# tex file sigma_13b

with open('build/sigma_13b.tex', 'w') as f:
  f.write(make_SI(sigma_13b,r'', figures=0))

# tex file rho

with open('build/rho.tex', 'w') as f:
  f.write(make_SI(rho,r'\per\meter\cubed', figures=2))

# tex file lambda_13b

with open('build/lambda_13b.tex', 'w') as f:
  f.write(make_SI(L_13b*1e-41,r'\meter', exp='1e+41', figures=2))

# tex file lambda_13b_parsec

with open('build/lambda_13b_parsec.tex', 'w') as f:
  f.write(make_SI(L_13b_parsec*1e-24,r'\text{pc}', exp='1e+24', figures=2))

#13c)--------------------------------------------------------------------------------------------------------------------------

L_13c = L_13b_parsec*3.26156

percent = 15*10**9/L_13c*100

# tex file lambda_13c

with open('build/lambda_13c.tex', 'w') as f:
  f.write(make_SI(L_13c*1e-25,r'\text{ly}', exp='1e+25', figures=2))

# tex file percent

with open('build/percent.tex', 'w') as f:
  f.write(make_SI(percent*1e+14,r'\percent', exp='1e-14', figures=2))

print(L_13c)