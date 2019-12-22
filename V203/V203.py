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

T = 373 # Temperatur für größer 1 bar

K_k1bar = const.convert_temperature(C_k1bar[:77], 'c', 'K')

K_g1bar = const.convert_temperature(C_g1bar[:77], 'c', 'K') 

pinP_k1bar = p_k1bar[:77] * 1e05
pinP_g1bar = p_g1bar[:77] * 1e05 

# Erstellung der Plots #################################################################################################################################################################

# Plot für kleiner 1 bar

park1b, covmk1b = np.polyfit(1/(const.R * K_k1bar), np.log(pinP_k1bar[:77]), deg=1, cov=True)
errk1b = np.sqrt(np.diag(covmk1b))

plt.plot(1/(const.R * K_k1bar), np.log(pinP_k1bar[:77]), 'kx', label='Messwerte')
x_plot = np.linspace(0.0003, 0.00042, 1000)
plt.plot(x_plot,park1b[1] +park1b[0]*x_plot, 'r-', label="Lineare Regression")
plt.xticks([3*1e-04, 3.25*1e-04, 3.5*1e-04, 3.75*1e-04,4 *1e-04],
           [3, 3.25, 3.5, 3.75, 4])
plt.legend(loc="best")
plt.xlabel(r'$1/(R \cdot T)*10^{-4}/(mol\, s^2/(kg \,m^2) $')
plt.ylabel(r' logarithmischer Dampfdruck $ln(p)$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotk1b.pdf')
plt.close()

park1b=unp.uarray(park1b, errk1b)

# Plot für größer 1 bar

parg1b, covmg1b = np.polyfit(K_g1bar, pinP_g1bar, deg=3, cov=True)
err = np.sqrt(np.diag(covmg1b))
print(parg1b)

VD = (const.R * K_g1bar)/(2*pinP_g1bar) - np.sqrt((const.R**2*K_g1bar**2)/(4*pinP_g1bar**2)-0.9/pinP_g1bar)
print(VD)

# Parameter werden ins tex-Format geschrieben ########################################################################################################################################################

with open('build/m_plotk1b.tex', 'w') as f: # es fehlen hier noch die Fehler
  f.write(make_SI(park1b[0]*1e-04,r'', figures=3)) # passt hier die Größenordnung?

with open('build/b_plotk1b.tex', 'w') as f: # es fehlen hier noch die Fehler
  f.write(make_SI(park1b[1],r'', figures=3))

# Tabellen ########################################################################################################################################################

# Tabelle für datak1b.tex -------------------------------------------------------------------------------------------------------------------------------------------------

p_k1bar1, p_k1bar2, p_k1bar3 = np.array_split(p_k1bar*1e+03, 3)
C_k1bar1, C_k1bar2, C_k1bar3 = np.array_split(C_k1bar, 3)

table_header = r'''
  \begin{tabular}{c c c c c c}
    \toprule
    \multicolumn{1}{c}{Druck} & \multicolumn{1}{c}{Temperatur} & \multicolumn{1}{c}{Druck} & \multicolumn{1}{c}{Temperatur} & \multicolumn{1}{c}{Druck} & \multicolumn{1}{c}{Temperatur}\\
    {in $\si{\milli\bar}$} & {in $\si{\celsius}$} &
    {in $\si{\milli\bar}$} & {in $\si{\celsius}$} &
    {in $\si{\milli\bar}$} & {in $\si{\celsius}$} \\

    \cmidrule(lr){1-2} \cmidrule(lr){3-4} \cmidrule(lr){5-6}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.0f} & {1:1.0f} & {2:1.0f} & {3:1.0f} & {4:1.0f} & {5:1.0f} \\'


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
    \multicolumn{1}{c}{Druck} & \multicolumn{1}{c}{Temperatur} & \multicolumn{1}{c}{Druck} & \multicolumn{1}{c}{Temperatur} & \multicolumn{1}{c}{Druck} & \multicolumn{1}{c}{Temperatur}\\
    {in $\si{\milli\bar}$} & {in $\si{\celsius}$} &
    {in $\si{\milli\bar}$} & {in $\si{\celsius}$} &
    {in $\si{\milli\bar}$} & {in $\si{\celsius}$} \\

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