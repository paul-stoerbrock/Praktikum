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

dnu15, rpm15 = np.genfromtxt('15degree.txt', unpack=True)
dnu30, rpm30 = np.genfromtxt('30degree.txt', unpack=True)
dnu45, rpm45 = np.genfromtxt('45degree.txt', unpack=True)
dnu, depth, strength = np.genfromtxt('data.txt', unpack=True)

#Konstanten ========================================================================================================================================================================================================================

#Sondenfrequenz:
nu_0 = 2*1e+6 #Hz

#Phantomflüssigkeit:
rho = 1.15  # g/ccm
c_L = 1800  # m/s
eta = 12    # mPa/s

#Dopplerprisma:
c_P = 2700  # m/s
l   = 30.7  # mm

#Formeln ========================================================================================================================================================================================================================

def alpha15(c_L, c_P):
  return np.pi/2-np.arcsin(np.sin(np.pi/6)*c_L/c_P)

def v15(nu_0, alpha15, dnu15, c_L):
  return 1/2*(np.abs(dnu15)*c_L)/(nu_0*np.cos(alpha15))

def alpha30(c_L, c_P):
  return np.pi/2-np.arcsin(np.sin(np.pi/3)*c_L/c_P)

def v30(nu_0, alpha30, dnu30, c_L):
  return 1/2*(np.abs(dnu30)*c_L)/(nu_0*np.cos(alpha30))

def alpha45(c_L, c_P):
  return np.pi/2-np.arcsin(np.sin(np.pi/2)*c_L/c_P)

def v45(nu_0, alpha45, dnu45, c_L):
  return 1/2*(np.abs(dnu45)*c_L)/(nu_0*np.cos(alpha45))

def alpha(c_L, c_P):
  return np.pi/2-np.arcsin(np.sin(np.pi/6)*c_L/c_P)

def v(nu_0, alpha, dnu, c_L):
  return 1/2*(np.abs(dnu)*c_L)/(nu_0*np.cos(alpha))


print(v15(nu_0, alpha15(c_L, c_P), dnu15, c_L))
print(v30(nu_0, alpha30(c_L, c_P), dnu30, c_L))
print(v45(nu_0, alpha45(c_L, c_P), dnu45, c_L))



#Graphen (WIP) ========================================================================================================================================================================================================================

par, covm = np.polyfit(np.abs(dnu15)/np.cos(alpha15(c_L, c_P)), v15(nu_0, alpha15(c_L, c_P), dnu15, c_L), deg=1, cov=True)
err = np.sqrt(np.diag(covm))

x_plot = np.linspace(250, 650, 1000)
plt.plot(x_plot ,par[0]*x_plot+par[1], 'b-', label="Lineare Regression")
plt.plot(np.abs(dnu15)/np.cos(alpha15(c_L, c_P)), v15(nu_0, alpha15(c_L, c_P), dnu15, c_L), 'ko', label="$Messwerte$")
plt.xlabel(r'$\frac{\Delta\nu}{\cos(\alpha)}$')
plt.ylabel(r'Strömungsgeschwindigkeit $v\;[\frac{m}{s}]$')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plot15.pdf')
plt.close()

m15 = ufloat(par[0], err[0])
b15 = ufloat(par[1], err[1])

# tex file for m15

with open('build/m15.tex', 'w') as f:
  f.write(make_SI(m15*1e+3,r'\milli\meter\per\hertz\per\second', figures=2))

# tex file for b15

with open('build/b15.tex', 'w') as f:
  f.write(make_SI(b15*1e+15,r'\femto\meter\per\second', figures=2))



par, covm = np.polyfit(np.abs(dnu30)/np.cos(alpha30(c_L, c_P)), v30(nu_0, alpha30(c_L, c_P), dnu30, c_L), deg=1, cov=True)
err = np.sqrt(np.diag(covm))

x_plot = np.linspace(200, 700, 1000)
plt.plot(x_plot ,par[0]*x_plot+par[1], 'b-', label="Lineare Regression")
plt.plot(np.abs(dnu30)/np.cos(alpha30(c_L, c_P)), v30(nu_0, alpha30(c_L, c_P), dnu30, c_L), 'ko', label="$Messwerte$")
plt.xlabel(r'$\frac{\Delta\nu}{\cos(\alpha)}$')
plt.ylabel(r'Strömungsgeschwindigkeit $v\;[\frac{m}{s}]$')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plot30.pdf')
plt.close()

m30 = ufloat(par[0], err[0])
b30 = ufloat(par[1], err[1])

# tex file for m30

with open('build/m30.tex', 'w') as f:
  f.write(make_SI(m30*1e+3,r'\milli\meter\per\hertz\per\second', figures=2))

# tex file for b30

with open('build/b30.tex', 'w') as f:
  f.write(make_SI(b30*1e+15,r'\femto\meter\per\second', figures=2))



par, covm = np.polyfit(np.abs(dnu45)/np.cos(alpha45(c_L, c_P)), v45(nu_0, alpha45(c_L, c_P), dnu45, c_L), deg=1, cov=True)
err = np.sqrt(np.diag(covm))

x_plot = np.linspace(350, 910, 1000)
plt.plot(x_plot ,par[0]*x_plot+par[1], 'b-', label="Lineare Regression")
plt.plot(np.abs(dnu45)/np.cos(alpha45(c_L, c_P)), v45(nu_0, alpha45(c_L, c_P), dnu45, c_L), 'ko', label="$Messwerte$")
plt.xlabel(r'$\frac{\Delta\nu}{\cos(\alpha)}$')
plt.ylabel(r'Strömungsgeschwindigkeit $v\;[\frac{m}{s}]$')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plot45.pdf')
plt.close()

m45 = ufloat(par[0], err[0])
b45 = ufloat(par[1], err[1])

# tex file for m45

with open('build/m45.tex', 'w') as f:
  f.write(make_SI(m45*1e+3,r'\milli\meter\per\hertz\per\second', figures=2))

# tex file for b45

with open('build/b45.tex', 'w') as f:
  f.write(make_SI(b45*1e+15,r'\femto\meter\per\second', figures=2))



#plt.plot(rpm45, dnu45, 'ko', label="$Messwerte$")
#plt.xlabel(r'Revolutions per Minute $[rpm]$')
#plt.ylabel(r'Frequenzverschiebung $\Delta\nu\;[nA]$')
#plt.legend(loc="best")
#plt.grid()
#plt.tight_layout
#plt.savefig('build/plot45.pdf')
#plt.close()
#
#
#
plt.plot(depth, v(nu_0, alpha(c_L, c_P), dnu, c_L), 'ko', label="$Messwerte$")
plt.axhline(y=max(v(nu_0, alpha(c_L, c_P), dnu, c_L)), color='r', linestyle= '--', label="Maximale Strömung")
plt.axvline(x=depth[7], color='b', linestyle= '--', label="Tiefe")
plt.xlabel(r'Tiefe $d\;[\mu s]$')
plt.ylabel(r'Strömungsgeschwindigkeit $v\;[\frac{m}{s}]$')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plotdepth.pdf')
plt.close()

plt.plot(strength, v(nu_0, alpha(c_L, c_P), dnu, c_L), 'ko', label="$Messwerte$")
plt.axhline(y=max(v(nu_0, alpha(c_L, c_P), dnu, c_L)), color='r', linestyle= '--', label="Maximale Strömung")
plt.axvline(x=strength[7], color='b', linestyle= '--', label="Stärkstes Signal")
plt.xlabel(r'Signalstärke $[\frac{V^2}{s}]$')
plt.ylabel(r'Strömungsgeschwindigkeit $v\;[\frac{m}{s}]$')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plotstrength.pdf')
plt.close()

#Tabellen ========================================================================================================================================================================================================================

#Messwerte ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

table_header = r'''
  \begin{tabular}{c c}
    \toprule
    {$\text{Frequenzverschiebung $\Delta\nu$} \; [\si{\hertz}]$} & {$\text{Revolutions per minute} \; [rpm]$} \\

    \cmidrule(lr{0,5em}){1-2}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.0f} & {1:1.0f} \\'

with open('build/15.tex', 'w') as g:
    g.write(table_header)
    for row in zip(dnu15, rpm15):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)



table_header = r'''
  \begin{tabular}{c c}
    \toprule
    {$\text{Frequenzverschiebung $\Delta\nu$} \; [\si{\hertz}]$} & {$\text{Revolutions per minute} \; [rpm]$} \\

    \cmidrule(lr{0,5em}){1-2}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.0f} & {1:1.0f} \\'

with open('build/30.tex', 'w') as g:
    g.write(table_header)
    for row in zip(dnu30, rpm30):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)



table_header = r'''
  \begin{tabular}{c c}
    \toprule
    {$\text{Frequenzverschiebung $\Delta\nu$} \; [\si{\hertz}]$} & {$\text{Revolutions per minute} \; [rpm]$} \\

    \cmidrule(lr{0,5em}){1-2}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.0f} & {1:1.0f} \\'

with open('build/45.tex', 'w') as g:
    g.write(table_header)
    for row in zip(dnu45, rpm45):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)



table_header = r'''
  \begin{tabular}{c c c}
    \toprule
    {$\text{Frequenzverschiebung $\Delta\nu$} \; [\si{\hertz}]$} & {$\text{Tiefe $d$} \; [\si{\micro\second}]$} & {$\text{Signalstärke} \; [\si{\volt\squared\per\second}]$} \\

    \cmidrule(lr{0,5em}){1-3}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.0f} & {1:1.1f} & {2:1.0f} \\'

with open('build/data.tex', 'w') as g:
    g.write(table_header)
    for row in zip(dnu, depth, strength):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)





#Geschwindigkeiten ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

table_header = r'''
  \begin{tabular}{c c}
    \toprule
    {$\text{Strömungsgeschwindigkeit $v$} \; [\si{\meter\per\second}]$} & {$\text{Revolutions per minute} \; [rpm]$} \\

    \cmidrule(lr{0,5em}){1-2}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.3f} & {1:1.0f} \\'

with open('build/v15.tex', 'w') as g:
    g.write(table_header)
    for row in zip( v15(nu_0, alpha15(c_L, c_P), dnu15, c_L), rpm15):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)



table_header = r'''
  \begin{tabular}{c c}
    \toprule
    {$\text{Strömungsgeschwindigkeit $v$} \; [\si{\meter\per\second}]$} & {$\text{Revolutions per minute} \; [rpm]$} \\

    \cmidrule(lr{0,5em}){1-2}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.3f} & {1:1.0f} \\'

with open('build/v30.tex', 'w') as g:
    g.write(table_header)
    for row in zip( v30(nu_0, alpha30(c_L, c_P), dnu30, c_L), rpm30):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)



table_header = r'''
  \begin{tabular}{c c}
    \toprule
    {$\text{Strömungsgeschwindigkeit $v$} \; [\si{\meter\per\second}]$} & {$\text{Revolutions per minute} \; [rpm]$} \\

    \cmidrule(lr{0,5em}){1-2}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.3f} & {1:1.0f} \\'

with open('build/v45.tex', 'w') as g:
    g.write(table_header)
    for row in zip( v45(nu_0, alpha45(c_L, c_P), dnu45, c_L), rpm45):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)



table_header = r'''
  \begin{tabular}{c c}
    \toprule
    {$\text{Strömungsgeschwindigkeit $v$} \; [\si{\meter\per\second}]$} & {$\text{Revolutions per minute} \; [rpm]$} \\

    \cmidrule(lr{0,5em}){1-2}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.3f} & {1:1.0f} \\'

with open('build/v.tex', 'w') as g:
    g.write(table_header)
    for row in zip( v(nu_0, alpha(c_L, c_P), dnu, c_L), rpm45):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)