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
        x = '{0:.{1:}f}'.format(num, figures)

    return r'\SI{{{}{}}}{{{}}}'.format(x, exp, unit)

# Definition der Variablen ########################################################################################################################################################

p_k1bar, C_k1bar = np.genfromtxt('datak1b.txt', unpack=True)
p_g1bar, C_g1bar = np.genfromtxt('datag1b.txt', unpack=True)

# Tabellen ########################################################################################################################################################

# Tabelle für datak1b.tex -------------------------------------------------------------------------------------------------------------------------------------------------

p_k1bar1, p_k1bar2, p_k1bar3 = np.array_split(p_k1bar, 3)
C_k1bar1, C_k1bar2, C_k1bar3 = np.array_split(C_k1bar, 3)

table_header = r'''
  \begin{tabular}{c c c c c c}
    \toprule
    {$Druck \:/\: \si{\bar}$} & {$Temperatur \:/\: \si{\celsius}$} &
    {$Druck \:/\: \si{\bar}$} & {$Temperatur \:/\: \si{\celsius}$} &
    {$Druck \:/\: \si{\bar}$} & {$Temperatur \:/\: \si{\celsius}$} \\

    \cmidrule(lr){1-2} \cmidrule(lr){3-4} \cmidrule(lr){5-6}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f} & {2:1.2f} & {3:1.2f} & {4:1.2f} & {5:1.2f} \\'


with open('build/table_k1bar.tex', 'w') as g:
    g.write(table_header)
    for row in zip(p_k1bar1, C_k1bar1, p_k1bar2, C_k1bar2, p_k1bar3, C_k1bar3):
        g.write(row_template.format(*row).replace('nan', ''))
        g.write('\n')
    g.write(table_footer)

# Tabelle für datag1b.tex -------------------------------------------------------------------------------------------------------------------------------------------------

p_g1bar1, p_g1bar2, p_g1bar3 = np.array_split(p_g1bar, 3)
C_g1bar1, C_g1bar2, C_g1bar3 = np.array_split(C_g1bar, 3)

table_header = r'''
  \begin{tabular}{c c c c c c}
    \toprule
    {$Druck \:/\: \si{\bar}$} & {$Temperatur \:/\: \si{\celsius}$} &
    {$Druck \:/\: \si{\bar}$} & {$Temperatur \:/\: \si{\celsius}$} &
    {$Druck \:/\: \si{\bar}$} & {$Temperatur \:/\: \si{\celsius}$} \\

    \cmidrule(lr){1-2} \cmidrule(lr{0,5em}){3-4} \cmidrule(lr{0,5em}){5-6}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.0f} & {2:1.1f} & {3:1.0f} & {4:1.1f} & {5:1.0f} \\'


with open('build/table_g1bar.tex', 'w') as g:
    g.write(table_header)
    for row in zip(p_g1bar1, C_g1bar1, p_g1bar2, C_g1bar2, p_g1bar3, C_g1bar3):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)