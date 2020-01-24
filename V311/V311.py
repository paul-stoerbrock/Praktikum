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

# Variablen ########################################################################################################################################################

Cu_IB, Cu_UHB, Cu_IQ, Cu_UHQ = np.genfromtxt('cu.txt', unpack=True)

Zn_IB, Zn_UHB, Zn_IQ, Zn_UHQ = np.genfromtxt('zn.txt', unpack=True)

Ag_IB, Ag_UHB, Ag_IQ, Ag_UHQ = np.genfromtxt('ag.txt', unpack=True)

hy_Iauf, hy_Bauf, hy_Iab, hy_Bab = np.genfromtxt('hysterese.txt', unpack=True)


# Plots ########################################################################################################################################################

#Plot von Hysteresekurve von der Spule
par1, covm1 = np.polyfit(hy_Iauf, hy_Bauf, deg=1, cov=True)
err1 = np.sqrt(np.diag(covm1))

par2, covm2 = np.polyfit(hy_Iab[1:9], hy_Bab[1:9], deg=1, cov=True)
err2 = np.sqrt(np.diag(covm2))


plt.plot(hy_Iauf, hy_Bauf,'kx', label='Messwerte')
plt.plot(hy_Iab, hy_Bab,'bx', label='Messwerte')
x_plot = np.linspace(0, 5, 1000)
plt.plot(x_plot, x_plot*par1[0]+par1[1], 'r-', label="Lineare Regression Hysteresekurve auf")
plt.plot(x_plot, x_plot*par2[0]+par2[1], 'g-', label="Lineare Regression Hysteresekurve ab")

plt.legend(loc="best")
plt.xlabel(r'Stromstärke $I\:/\:V$')
plt.ylabel(r'magnetische Feldstärke $B\:/\:T$')
plt.grid()
plt.tight_layout
plt.savefig('build/plothy.pdf')
plt.close()

parhyauf=unp.uarray(par1, err1)
parhyab=unp.uarray(par2, err2)

#Plot von Kupfer mit konstantem Querstrom
B= (parhyauf[0].n + parhyab[0].n)/2*Cu_IB

par, covm = np.polyfit(B, Cu_UHB, deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(B, Cu_UHB,'kx', label='Messwerte')
x_plot = np.linspace(0, 1.3, 1000)
plt.plot(x_plot, x_plot*par[0]+par[1], 'r-', label="Lineare Regression")
plt.yticks([0, 5*1e-06, 10*1e-06, 15*1e-06, 20*1e-06],
           [0, 5, 10, 15, 20])
plt.legend(loc="best")
plt.xlabel(r'magnetische Feldstärke $B\:/\:T$')
plt.ylabel(r'Hallspannung $U_H\:/\:\mu V$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotCu_IB.pdf')
plt.close()

parCu_IB=unp.uarray(par, err)

#Plot von Kupfer mit konstantem B-Feld
par, covm = np.polyfit(Cu_IQ, Cu_UHQ, deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(Cu_IQ, Cu_UHQ,'kx', label='Messwerte')
x_plot = np.linspace(0, 10, 1000)
plt.plot(x_plot, x_plot*par[0]+par[1], 'r-', label="Lineare Regression")
plt.yticks([0, 5*1e-06, 10*1e-06, 15*1e-06, 20*1e-06],
           [0, 5, 10, 15, 20])
plt.legend(loc="best")
plt.xlabel(r'Quellstrom $I_Q\:/\:A$')
plt.ylabel(r'Hallspannung $U_H\:/\:\mu V$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotCu_IQ.pdf')
plt.close()

parCu_IQ=unp.uarray(par, err)

#Plot von Silber mit konstantem Querstrom
B= (parhyauf[0].n + parhyab[0].n)/2*Ag_IB

par, covm = np.polyfit(B, Ag_UHB, deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(B, Ag_UHB,'kx', label='Messwerte')
x_plot = np.linspace(0, 1.3, 1000)
plt.plot(x_plot, x_plot*par[0]+par[1], 'r-', label="Lineare Regression")
#plt.yticks([0, 5*1e-06, 10*1e-06, 15*1e-06, 20*1e-06],
#           [0, 5, 10, 15, 20])
plt.legend(loc="best")
plt.xlabel(r'magnetische Feldstärke $B\:/\:T$')
plt.ylabel(r'Hallspannung $U_H\:/\:\mu V$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotAg_IB.pdf')
plt.close()

parAg_IB=unp.uarray(par, err)

#Plot von Silber mit konstantem B-Feld
par, covm = np.polyfit(Ag_IQ, Ag_UHQ, deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(Ag_IQ, Ag_UHQ,'kx', label='Messwerte')
x_plot = np.linspace(0, 10, 1000)
plt.plot(x_plot, x_plot*par[0]+par[1], 'r-', label="Lineare Regression")
#plt.yticks([0, 5*1e-06, 10*1e-06, 15*1e-06, 20*1e-06],
#           [0, 5, 10, 15, 20])
plt.legend(loc="best")
plt.xlabel(r'Quellstrom $I_Q\:/\:A$')
plt.ylabel(r'Hallspannung $U_H\:/\:\mu V$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotAg_IQ.pdf')
plt.close()

parAg_IQ=unp.uarray(par, err)


# Tabellen ########################################################################################################################################################

# Kupfer ===============================================================================================================================================================

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    \multicolumn{2}{c}{$\text{Konstanter Querstrom =} \; \SI{10}{\ampere}$} & \multicolumn{2}{c}{$\text{Konstantes B-Feld =} \; \SI{5}{\ampere}$}\\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4}
    \multicolumn{1}{c}{$\text{Strom} \; I_B \:/\: \si{\ampere}$} & \multicolumn{1}{c}{$\text{Spannung} \; U_H \:/\: \si{\milli\volt} $} & 
    \multicolumn{1}{c}{$\text{Strom} \; I_Q \:/\: \si{\ampere}$} & \multicolumn{1}{c}{$\text{Spannung} \; U_H \:/\: \si{\milli\volt}$}\\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4} 

'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.2f} & {2:1.0f} & {3:1.2f}\\'

with open('build/Cu_table.tex', 'w') as g:
    g.write(table_header)
    for row in zip(Cu_IB, Cu_UHB*1e+04, Cu_IQ, Cu_UHQ*1e+04):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

# Zink ===============================================================================================================================================================

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    \multicolumn{2}{c}{$\text{Konstanter Querstrom =} \; \SI{8}{\ampere}$} & \multicolumn{2}{c}{$\text{Konstantes B-Feld =} \; \SI{5}{\ampere}$}\\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4}
    \multicolumn{1}{c}{$\text{Strom} \; I_B \:/\: \si{\ampere}$} & \multicolumn{1}{c}{$\text{Spannung} \; U_H \:/\: \si{\milli\volt} $} & 
    \multicolumn{1}{c}{$\text{Strom} \; I_Q \:/\: \si{\ampere}$} & \multicolumn{1}{c}{$\text{Spannung} \; U_H \:/\: \si{\milli\volt}$}\\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4} 

'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.2f} & {2:1.1f} & {3:1.2f}\\'

with open('build/Zn_table.tex', 'w') as g:
    g.write(table_header)
    for row in zip(Zn_IB, Zn_UHB*1e+04, Zn_IQ, Zn_UHQ*1e+04):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

# Silber ===============================================================================================================================================================

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    \multicolumn{2}{c}{$\text{Konstanter Querstrom =} \; \SI{10}{\ampere}$} & \multicolumn{2}{c}{$\text{Konstantes B-Feld =} \; \SI{5}{\ampere}$}\\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4}
    \multicolumn{1}{c}{$\text{Strom} \; I_B \:/\: \si{\ampere}$} & \multicolumn{1}{c}{$\text{Spannung} \; U_H \:/\: \si{\milli\volt} $} & 
    \multicolumn{1}{c}{$\text{Strom} \; I_Q \:/\: \si{\ampere}$} & \multicolumn{1}{c}{$\text{Spannung} \; U_H \:/\: \si{\milli\volt}$}\\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4} 

'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.2f} & {2:1.0f} & {3:1.2f}\\'

with open('build/Ag_table.tex', 'w') as g:
    g.write(table_header)
    for row in zip(Ag_IB, Ag_UHB*1e+04, Ag_IQ, Ag_UHQ*1e+04):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

# Hysterese ===============================================================================================================================================================

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    \multicolumn{2}{c}{$\text{Aufsteigend von} \; \SI{0}{\ampere}$} & \multicolumn{2}{c}{$\text{Absteigend von} \; \SI{4.5}{\ampere}$}\\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4}
    \multicolumn{1}{c}{$\text{Strom} \; I_{auf} \:/\: \si{\ampere}$} & \multicolumn{1}{c}{$\text{Magnetfeld} \; B_{auf} \:/\: \si{\milli\tesla} $} & 
    \multicolumn{1}{c}{$\text{Strom} \; I_{ab} \:/\: \si{\ampere}$} & \multicolumn{1}{c}{$\text{Magnetfeld} \; B_{ab} \:/\: \si{\milli\tesla}$}\\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4} 

'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.1f} & {2:1.1f} & {3:1.1f}\\'

with open('build/Hysterese_table.tex', 'w') as g:
    g.write(table_header)
    for row in zip(hy_Iauf, hy_Bauf*1e+03, hy_Iab, hy_Bab*1e+03):
        g.write(row_template.format(*row).replace('nan',' ') )
        g.write('\n')
    g.write(table_footer)