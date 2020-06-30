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

lit = const.h/const.e





#Graphen ========================================================================================================================================================================================================================

par, covm = np.polyfit(-U_Gelb[12:26], np.sqrt(I_Gelb[12:26]), deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(-U_Gelb, np.sqrt(I_Gelb), 'kx', label='Messwerte')
x_plot = np.linspace(-7, 0, 1000)
plt.plot(x_plot ,par[0]*x_plot+par[1], 'b-', label="Lineare Regression")

U_g_gelb = -par[1]/par[0]
plt.plot(U_g_gelb, 0, 'ko', label="$U_g$")
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
U_g_gelb_err = -b_gelb/m_gelb

# tex file for m_gelb 

with open('build/m_gelb.tex', 'w') as f:
  f.write(make_SI(m_gelb,r'\nano\ampere\tothe{1/2}\per\volt', figures=2))

# tex file for b_gelb 

with open('build/b_gelb.tex', 'w') as f:
  f.write(make_SI(b_gelb,r'\nano\ampere', figures=2))

# tex file for U_g_gelb 

with open('build/U_g_gelb.tex', 'w') as f:
  f.write(make_SI(U_g_gelb,r'\volt', figures=2))

# tex file for U_g_gelb_err 

with open('build/U_g_gelb_err.tex', 'w') as f:
  f.write(make_SI(U_g_gelb_err,r'\volt', figures=2))





par, covm = np.polyfit(-U_Gruen, np.sqrt(I_Gruen), deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(-U_Gruen, np.sqrt(I_Gruen), 'kx', label='Messwerte')
x_plot = np.linspace(-2, 0.25, 1000)
plt.plot(x_plot ,par[0]*x_plot+par[1], 'b-', label="Lineare Regression")

#0=par[0]*x_plot+par[1]
U_g_gruen = -par[1]/par[0]
plt.plot(U_g_gruen, 0, 'ko', label="$U_g$")
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
  f.write(make_SI(m_gruen,r'\nano\ampere\tothe{1/2}\per\volt', figures=2))

# tex file for b_gruen 

with open('build/b_gruen.tex', 'w') as f:
  f.write(make_SI(b_gruen,r'\nano\ampere', figures=2))

# tex file for U_g_gruen 

with open('build/U_g_gruen.tex', 'w') as f:
  f.write(make_SI(U_g_gruen,r'\volt', figures=2))





par, covm = np.polyfit(-U_Rot, np.sqrt(I_Rot), deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(-U_Rot, np.sqrt(I_Rot), 'kx', label='Messwerte')
x_plot = np.linspace(-2, 1.5, 1000)
plt.plot(x_plot ,par[0]*x_plot+par[1], 'b-', label="Lineare Regression")

U_g_rot = -par[1]/par[0]
plt.plot(U_g_rot, 0, 'ko', label="$U_g$")
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
  f.write(make_SI(m_rot,r'\nano\ampere\tothe{1/2}\per\volt', figures=2))

# tex file for b_rot 

with open('build/b_rot.tex', 'w') as f:
  f.write(make_SI(b_rot,r'\nano\ampere', figures=2))

# tex file for U_g_rot 

with open('build/U_g_rot.tex', 'w') as f:
  f.write(make_SI(U_g_rot,r'\volt', figures=2))



f_gelb = const.c/578e-09
f_gruen = const.c/546e-09
f_rot = const.c/671e-09

# tex file for f_gelb

with open('build/f_gelb.tex', 'w') as f:
  f.write(make_SI(f_gelb*1e-15,r'\femto\hertz', figures=2))

# tex file for f_gruen

with open('build/f_gruen.tex', 'w') as f:
  f.write(make_SI(f_gruen*1e-15,r'\femto\hertz', figures=2))

# tex file for f_rot

with open('build/f_rot.tex', 'w') as f:
  f.write(make_SI(f_rot*1e-15,r'\femto\hertz', figures=2))

U_g_array = np.array([U_g_gelb, U_g_gruen, U_g_rot])
f_array = np.array([f_gelb, f_gruen, f_rot])

par, covm = np.polyfit(f_array, U_g_array, deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(f_array, U_g_array, 'kx', label='Messwerte')
x_plot = np.linspace(0e+14, 6e+14, 1000)
plt.plot(x_plot ,par[0]*x_plot+par[1], 'b-', label="Lineare Regression")

plt.xlabel(r'Frequenz $\nu\;[Hz]$')
plt.ylabel(r'Spannung $U_g\;[V]$')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plotLambda.pdf')
plt.close()

m_f = ufloat(par[0], err[0])
b_f = ufloat(par[1], err[1])
fehler = np.abs(m_f.n-lit)/lit

# tex file for m_f 

with open('build/m_f.tex', 'w') as f:
  f.write(make_SI(m_f*1e+14,r'\volt\per\hertz', figures=2))

# tex file for b_f 

with open('build/b_f.tex', 'w') as f:
  f.write(make_SI(b_f,r'\volt', figures=2))

# tex file for fehler 

with open('build/fehler.tex', 'w') as f:
  f.write(make_SI(fehler*100,r'\percent', figures=2))

print(fehler)
print(b_f)
print(m_f)
print(b_f)
print(m_f)

U_g_array_L = np.array([U_g_gelb, U_g_gruen])
f_array_L = np.array([f_gelb, f_gruen])


plt.plot(f_array_L, U_g_array_L, 'kx', label='Messwerte')
m_L = (U_g_gruen-U_g_gelb)/(f_gruen-f_gelb)
b_L =U_g_gelb-m_L*f_gelb
f_array_L = np.linspace(0e+14, 6e+14, 1000)
plt.plot(f_array_L, m_L*f_array_L+b_L, 'b-', label="Gerade")

plt.xlabel(r'Frequenz $\nu\;[Hz]$')
plt.ylabel(r'Spannung $U_g\;[V]$')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plotLambdaL.pdf')
plt.close()

fehler_L = np.abs(m_L-lit)/lit

# tex file for m_L 

with open('build/m_L.tex', 'w') as f:
  f.write(make_SI(m_L*1e+14,r'\volt\per\hertz', figures=2))

# tex file for b_L 

with open('build/b_L.tex', 'w') as f:
  f.write(make_SI(b_L,r'\volt', figures=2))

# tex file for fehler_L 

with open('build/fehler_L.tex', 'w') as f:
  f.write(make_SI(fehler_L*100,r'\percent', figures=2))

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