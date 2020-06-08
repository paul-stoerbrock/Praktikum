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

h, I = np.genfromtxt('data.txt', unpack=True)








def I_h(x, a, b, c):
  return a*x+b*x**2+c*x**3

popT, pcov = curve_fit(
    I_h,
    h,
    I,
    sigma=None,
    absolute_sigma=True,
    p0=[100, 28, -30]
    )
err_var = np.sqrt(np.diag(pcov))

plt.plot(h, I, 'kx', label='Messwerte')
x_plot = np.linspace(0, 10000, 100000)
plt.plot(x_plot, I_h(x_plot,*popT), linestyle='-', label='Ausgleichskurve')

plt.xlabel(r'Höhe h [m]')
plt.ylabel(r'Intensität I [$cm^{-3}s^{-1}$]')
plt.axvline(x=1000, color='g', linestyle= ':', label='$I_{min}$')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plot18.pdf')
plt.close()

I_min = np.min(I_h(x_plot,*popT))
print(I_min)
h_min = h[1]
print(h_min)

a = popT[0]
b = popT[1] 
c = popT[2]

# tex file for a

with open('build/a.tex', 'w') as f:
  f.write(make_SI(a*1e+03,r'\per\centi\meter\squared\per\second', exp='e-03', figures=4))

# tex file for b

with open('build/b.tex', 'w') as f:
  f.write(make_SI(b*1e+07,r'\per\centi\meter\per\second', exp='e-07', figures=4))

# tex file for c

with open('build/c.tex', 'w') as f:
  f.write(make_SI(c*1e+11,r'\per\cubed\per\second', exp='e-11', figures=4))

# tex file h_min

with open('build/h_min.tex', 'w') as f:
  f.write(make_SI(h_min,r'\meter', figures=0))





table_header = r'''
  \begin{tabular}{c c}
    \toprule
    {$\text{Höhe h} \; [\si{\meter}]$} & {$\text{Intensität I} \; [\si{\per\centi\meter\cubed\per\second}]$} \\

    \cmidrule(lr{0,5em}){1-2}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.0f} & {1:1.2f}\\'

with open('build/table18.tex', 'w') as g:
    g.write(table_header)
    for row in zip(h,I):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)