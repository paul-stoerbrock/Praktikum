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

U_GM, N_GM = np.genfromtxt('Kennlinie.dat', unpack = True)
U_strom, I_strom = np.genfromtxt('Zaehlrohrstrom.dat', unpack = True)

N_1 = 96041/2 #per second
N_2 = 76518/2 #per second
N_12 = 158479/2 #per second

T=(N_1+N_2-N_12)/(2*N_1*N_2) #micro second

N_GM_error = unp.uarray(N_GM, np.sqrt(N_GM))
I_bar_err = unp.uarray(I_strom, 0.05) #in micro Amp

N_array =np.array([N_GM[3], N_GM[8], N_GM[13], N_GM[18], N_GM[23], N_GM[28], N_GM[33], N_GM[38]])
N_array_err = unp.uarray(N_array, np.sqrt(N_array))

Z= (I_bar_err*1e-06)/(const.e*N_array_err)

#Kennlinie (Plateau)---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
par, covm = np.polyfit(U_GM[3:35], N_GM[3:35], deg=1, cov=True)
err = np.sqrt(np.diag(covm))

yerr=np.sqrt(N_GM)

plt.errorbar(U_GM, N_GM, yerr, fmt='kx')
x_plot = np.linspace(350, 670, 1000)
plt.plot(x_plot ,par[0]*x_plot+par[1], 'b-', label="Lineare Regression des Plateaus")
plt.xlabel(r'Spannung$\;U\;[V]$')
plt.ylabel(r'Impulse$\;N\;$[s$^{-1}$]')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plotGM.pdf')
plt.close()

m_GM = ufloat(par[0], err[0])
b_GM = ufloat(par[1], err[1])

#Zahl der freigesetzten Ladung---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
par, covm = np.polyfit(noms(I_bar_err*1e-06), noms(Z), deg=1, cov=True)
err = np.sqrt(np.diag(covm))

yerr=0.05*1e-06
xerr=stds(Z)

plt.plot(noms(I_bar_err*1e-06), noms(Z))
plt.errorbar(noms(I_bar_err*1e-06), noms(Z), xerr, yerr, fmt='kx')
x_plot = np.linspace(0.125*1e-06, 2*1e-06, 1000)
plt.plot(x_plot ,par[0]*x_plot+par[1], 'b-', label="Lineare Regression")
plt.xlabel(r'Strom$\;I\;$[Amp]')
plt.ylabel(r'Z')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plotZ.pdf')
plt.close()

print(stds(Z))

## tex file for theta_bragg_lit
#
#with open('build/theta_bragg_lit.tex', 'w') as f:
#  f.write(make_SI(theta_bragg_lit,r'°', figures=1))

#Geiger-Müller Kennlinie-------------------------------------------------------------------------------------------------------------------------------------------------

table_header = r'''
\begin{longtable}{S[table-format=2.1] S[table-format=3.0] S[table-format=4.1] S[table-format=2.0] S[table-format=3.1] S[table-format=4.0]}
    \caption{Messwerte der Geiger-Müller Kennlinie}\\
    \toprule
    {$\text{Spannung $U$} \; [\si{\volt}]$} & {$\text{Impulse $N$} \; [\si{\per\second}]$} \\ 
    \cmidrule(lr{0.5em}){1-2}
'''
table_footer = r''' 
    \bottomrule
    \label{tab:1}
\end{longtable}
'''

row_template = r'     {0:2.0f} & {1:3.0f} \\'

with open('build/GM.tex', 'w') as g:
    g.write(table_header)
    for row in zip(U_GM, N_GM):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

#Geiger-Müller Kennlinie-------------------------------------------------------------------------------------------------------------------------------------------------

#U_GM1, U_GM2, U_GM3 = np.array_split(U_GM, 3)
#N_GM1, N_GM2, N_GM3 = np.array_split(N_GM, 3)
#
#table_header = r'''
#  \begin{tabular}{c c c c c c}
#    \toprule
#    {$\text{Spannung $U$} \; [\si{\volt}]$} & {$\text{Impulse $N$} \; [\si{\per\second}]$} &
#    {$\text{Spannung $U$} \; [\si{\volt}]$} & {$\text{Impulse $N$} \; [\si{\per\second}]$} &
#    {$\text{Spannung $U$} \; [\si{\volt}]$} & {$\text{Impulse $N$} \; [\si{\per\second}]$} \\
#
#    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4} \cmidrule(lr{0,5em}){5-6}
#'''
#table_footer = r'''    \bottomrule
#  \end{tabular}
#'''
#row_template = r'     {0:1.2f} & {1:1.2f} & {2:1.2f} & {3:1.2f} & {4:1.2f} & {5:1.2f}  \\'
#
#with open('build/GM.tex', 'w') as g:
#    g.write(table_header)
#    for row in zip(U_GM1, N_GM1, U_GM2, N_GM2, U_GM3, N_GM3):
#        g.write(row_template.format(*row))
#        g.write('\n')
#    g.write(table_footer)

#Geiger-Müller Kennlinie-------------------------------------------------------------------------------------------------------------------------------------------------

table_header = r'''
  \begin{tabular}{c c c c c c}
    \toprule
    {$\text{Spannung $U$} \; [\si{\volt}]$} & {$\text{Strom $I$} \; [\si{\micro\volt}]$} \\

    \cmidrule(lr{0,5em}){1-2}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f} \\'

with open('build/strom.tex', 'w') as g:
    g.write(table_header)
    for row in zip(U_strom, I_strom):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)