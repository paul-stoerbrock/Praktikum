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


def I_f( phi, A0, b):
    return  A0**2 *np.sinc(b*np.sin(phi)/(532*1e-09))**2

#Daten ########################################################################

b, I = np.genfromtxt('data.txt', unpack=True)

b= b*1e-03
I_ohne_dunkel = I*1e-06-0.18*1e-09

L = ufloat(39.3, 0.1)*1e-02

phi = (b-27.7*1e-03)/(39.3*1e-02) 

# Plots



par, cov = curve_fit(
    I_f,
    noms(phi),
    noms(I_ohne_dunkel),
    sigma=None,
    absolute_sigma=True,
    p0=[ 0.008, 0.1*1e-03]
    )

err=np.sqrt(np.diag(cov))

plt.plot(phi, I_ohne_dunkel ,'kx' , label='Messwerte')
x_plot = np.linspace(-0.03, 0.06, 1000)
plt.plot(x_plot, I_f(x_plot, *par), 'r-', label='Fitkurve')
plt.legend(loc="best")
plt.yticks(
    [0, 1e-05, 2e-05, 3e-05, 4e-05, 5e-05, 6e-05, 7e-05],
    [0, 10, 20, 30, 40, 50, 60, 70]
)

par = unp.uarray(par,err)
print(par)
plt.xlabel(r' $\varphi \:/\:rad$')
plt.ylabel(r' $I\:/\:\mu A$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot.pdf')
plt.close()







# tex file for 

with open('build/A.tex', 'w') as f:
  f.write(make_SI(par[0],r'', figures=1))

with open('build/b.tex', 'w') as f:
  f.write(make_SI(par[1]*1e03,r'\milli\meter', figures=3))

with open('build/relerr_b.tex', 'w') as f:
  f.write(make_SI(rel_err(par[1], 0.1*1e-03) ,r'\percent', figures=2))



# Tabellen #######################################################################

I1, I2 = np.array_split(I, 2)
b1, b2 = np.array_split(b, 2)

table_header = r'''
  \begin{longtable}[H]{S[table-format=2.1] S[table-format=1.2] S[table-format=2.1] S[table-format=2.2]  }
    \caption{
        In der Tabelle sind die Messwerte aus dem Versuch mit dem Einzelspalt notiert.
        Gemessen wird die Verschiebung auf dem Verschiebereiter $b$ gegen die vom Photoelement
        gemessene Stromstärke I. 
    }\\ 
    \toprule
    \multicolumn{1}{c}{ $b\:/\:mm$ } & \multicolumn{1}{c}{$I\:/\:\mu A$ }   &
    \multicolumn{1}{c}{ $b\:/\:mm$ } & \multicolumn{1}{c}{$I\:/\:\mu A$ }   \\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4}  

'''
table_footer = r'''    \bottomrule
  \label{tab:1}
  \end{longtable}
'''
row_template = r'     {0:1.1f} & {1:1.2f} & {2:1.1f} & {3:1.2f}   \\'

with open('build/table.tex', 'w') as g:
    g.write(table_header)
    for row in zip( b1*1e03, I1, b2*1e03, I2):
        g.write(row_template.format(*row).replace('nan',' ').replace('.',',') )
        g.write('\n')
    g.write(table_footer)