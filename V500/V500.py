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

U_Rot, I_Rot = np.genfromtxt('Rot.txt', unpack=True)
U_Gelb, I_Gelb = np.genfromtxt('Gelb.txt', unpack=True)
U_Gruen, I_Gruen = np.genfromtxt('Gruen.txt', unpack=True)






#Graphen ========================================================================================================================================================================================================================

par, covm = np.polyfit(U_Gelb[12:26], np.sqrt(I_Gelb[12:26]), deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(U_Gelb, np.sqrt(I_Gelb), 'kx', label='Messwerte')
x_plot = np.linspace(0, 7, 1000)
plt.plot(x_plot ,par[0]*x_plot+par[1], 'b-', label="Lineare Regression")

U_g_gelb = -par[1]/par[0]
plt.plot(U_g_gelb, 0, 'ko', label="U_g")
plt.axhline(y=0, color='r', linestyle= '--', label="U-Achse")
plt.xlabel(r'Spannung $U\;[V]$')
plt.ylabel(r'Strom $\sqrt{I}\;[nA]$')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plotGelb.pdf')
plt.close()

m_gelb = ufloat(par[0],err[0])
b_gelb = ufloat(par[1],err[1])

# tex file for m_gelb 

with open('build/m_gelb.tex', 'w') as f:
  f.write(make_SI(m_gelb,r'\nano\ampere\tothe{1/2}\per\volt', figures=1))

# tex file for b_gelb 

with open('build/b_gelb.tex', 'w') as f:
  f.write(make_SI(b_gelb,r'\nano\ampere', figures=1))

# tex file for U_g_gelb 

with open('build/U_g_gelb.tex', 'w') as f:
  f.write(make_SI(U_g_gelb,r'\volt', figures=1))





par, covm = np.polyfit(U_Gruen, np.sqrt(I_Gruen), deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(U_Gruen, np.sqrt(I_Gruen), 'kx', label='Messwerte')
x_plot = np.linspace(-0.25, 2, 1000)
plt.plot(x_plot ,par[0]*x_plot+par[1], 'b-', label="Lineare Regression")

#0=par[0]*x_plot+par[1]
U_g_gruen = -par[1]/par[0]
plt.plot(U_g_gruen, 0, 'ko', label="U_g")
plt.axhline(y=0, color='r', linestyle= '--', label="U-Achse")
plt.xlabel(r'Spannung $U\;[V]$')
plt.ylabel(r'Strom $\sqrt{I}\;[nA]$')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plotGruen.pdf')
plt.close()

m_gruen = ufloat(par[0],err[0])
b_gruen = ufloat(par[1],err[1])

# tex file for m_gruen 

with open('build/m_gruen.tex', 'w') as f:
  f.write(make_SI(m_gruen,r'\nano\ampere\tothe{1/2}\per\volt', figures=1))

# tex file for b_gruen 

with open('build/b_gruen.tex', 'w') as f:
  f.write(make_SI(b_gruen,r'\nano\ampere', figures=1))

# tex file for U_g_gruen 

with open('build/U_g_gruen.tex', 'w') as f:
  f.write(make_SI(U_g_gruen,r'\volt', figures=1))





par, covm = np.polyfit(U_Rot, np.sqrt(I_Rot), deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(U_Rot, np.sqrt(I_Rot), 'kx', label='Messwerte')
x_plot = np.linspace(-1.5, 2, 1000)
plt.plot(x_plot ,par[0]*x_plot+par[1], 'b-', label="Lineare Regression")

U_g_rot = -par[1]/par[0]
plt.plot(U_g_rot, 0, 'ko', label="U_g")
plt.axhline(y=0, color='r', linestyle= '--', label="U-Achse")
plt.xlabel(r'Spannung $U\;[V]$')
plt.ylabel(r'Strom $\sqrt{I}\;[nA]$')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plotRot.pdf')
plt.close()

m_rot = ufloat(par[0],err[0])
b_rot = ufloat(par[1],err[1])

# tex file for m_rot 

with open('build/m_rot.tex', 'w') as f:
  f.write(make_SI(m_rot,r'\nano\ampere\tothe{1/2}\per\volt', figures=1))

# tex file for b_rot 

with open('build/b_rot.tex', 'w') as f:
  f.write(make_SI(b_rot,r'\nano\ampere', figures=1))

# tex file for U_g_rot 

with open('build/U_g_rot.tex', 'w') as f:
  f.write(make_SI(U_g_rot,r'\volt', figures=1))



#Tabellen ========================================================================================================================================================================================================================

table_header = r'''
  \begin{tabular}{c c}
    \toprule
    {$\text{Spannung $U$} \; [\si{\volt}]$} & {$\text{Strom $I$} \; [\si{\nano\ampere}]$} \\

    \cmidrule(lr{0,5em}){1-2}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f} \\'

with open('build/Gelb.tex', 'w') as g:
    g.write(table_header)
    for row in zip(U_Gelb, I_Gelb):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

table_header = r'''
  \begin{tabular}{c c}
    \toprule
    {$\text{Spannung $U$} \; [\si{\volt}]$} & {$\text{Strom $I$} \; [\si{\nano\ampere}]$} \\

    \cmidrule(lr{0,5em}){1-2}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f} \\'

with open('build/Rot.tex', 'w') as g:
    g.write(table_header)
    for row in zip(U_Rot, I_Rot):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

table_header = r'''
  \begin{tabular}{c c}
    \toprule
    {$\text{Spannung $U$} \; [\si{\volt}]$} & {$\text{Strom $I$} \; [\si{\nano\ampere}]$} \\

    \cmidrule(lr{0,5em}){1-2}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f} \\'

with open('build/Gruen.tex', 'w') as g:
    g.write(table_header)
    for row in zip(U_Gruen, I_Gruen):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

## tex file for 
#
#with open('build/.tex', 'w') as f:
#  f.write(make_SI(,r'', figures=1))