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

hy_Iauf, hy_Bauf, hy_Iab, hy_Bab = np.genfromtxt('hysterese.txt', unpack=True) # Es müssen für hy_Iab und hy_Bab die 'nan' entfernt werden!

# Plots ########################################################################################################################################################



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
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)