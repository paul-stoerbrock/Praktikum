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

theta_bragg, N_bragg = np.genfromtxt('Bragg.dat', unpack = True)
theta_brom, N_brom = np.genfromtxt('Brom.dat', unpack = True)
theta_cu, N_cu = np.genfromtxt('Emissionsspektrum.dat', unpack = True)
theta_gallium, N_gallium = np.genfromtxt('Gallium.dat', unpack = True)
theta_rub, N_rub = np.genfromtxt('Rubidium.dat', unpack = True)
theta_stron, N_stron = np.genfromtxt('Strontium.dat', unpack = True)
theta_zink, N_zink = np.genfromtxt('Zink.dat', unpack = True)
theta_zirk, N_zirk = np.genfromtxt('Zirkonium.dat', unpack = True)

#Globale Variablen---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

d_LiF = 201.4*1e-12 #piko meter
E_Kedge_cu = 8987.9615 #eV
R_infty = 13.6 #eV

#Ordnungszahlen---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

z_cu = 29

#Emissionsspektrum Kupfer---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

par, covm = np.polyfit(theta_cu, N_cu, deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(theta_cu, N_cu, 'kx', label='Messwerte')
x_plot = np.linspace(8, 30, 1000)

plt.xlabel(r'Winkel $\theta\;[\degree]$')
plt.ylabel(r'Impuls [s$^{-1}$]')
ymax_cu = np.max(N_cu)
xmax_cu = theta_cu[145]
ylmax_cu = N_cu[122]
xlmax_cu = theta_cu[122]
plt.axvline(x=xmax_cu, color='r', linestyle= ':', label='$K_{\\alpha}$-Linie')
plt.axvline(x=xlmax_cu, color='b', linestyle= ':', label='$K_{\\beta}$-Linie')
plt.annotate('$K_{\\alpha}$', xy=(xmax_cu, 2500), xycoords='data', xytext=(xmax_cu+1, 2500), textcoords='data', arrowprops=dict(arrowstyle='->', facecolor='grey'), horizontalalignment='left')
plt.annotate('$K_{\\beta}$', xy=(xlmax_cu, 2500), xycoords='data', xytext=(xlmax_cu+1, 2500), textcoords='data', arrowprops=dict(arrowstyle='->', facecolor='grey'), horizontalalignment='left')
plt.annotate('Bremsstrahlung', xy=(xmax_cu/2, 1100), xytext=(xmax_cu/2, 1100), verticalalignment='top', horizontalalignment='left')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plotCu.pdf')
plt.close()

xmax_cu_float = ufloat(xmax_cu, 0.1)
xlmax_cu_float = ufloat(xlmax_cu, 0.1)
xmax_cu_bogen = unp.deg2rad(xmax_cu_float)
xlmax_cu_bogen = unp.deg2rad(xlmax_cu_float)
Breite_lKa_cu = ufloat(20, 0.1)
Breite_rKa_cu = ufloat(20.6, 0.1)
BrBogen_lKa_cu = unp.deg2rad(Breite_lKa_cu)
BrBogen_rKa_cu = unp.deg2rad(Breite_rKa_cu)

HWB_cu = BrBogen_rKa_cu-BrBogen_lKa_cu

print(xmax_cu_F)


L_cu_a = 2*d_LiF*unp.sin(xmax_cu_F)
L_cu_b = 2*d_LiF*unp.sin(xlmax_cu_F)
L_HWB_cu = 2*d_LiF*unp.sin(np.deg2rad(HWB_cu))

print(L_cu_a)

E_cu_a = (const.h*const.c)/(L_cu_a*const.e)
E_cu_b = (const.h*const.c)/(L_cu_b*const.e)
E_HWB_cu = (const.h*const.c)/(L_HWB_cu*const.e)

print(E_cu_a)
print(E_HWB_cu)

sigma1_cu = 29-unp.sqrt(E_Kedge_cu/R_infty)
sigma2_cu = 29-unp.sqrt((4*(E_cu_a-E_Kedge_cu))/R_infty)
sigma3_cu = 29-unp.sqrt((9*(E_cu_a-E_Kedge_cu))/R_infty)

#Tabelle Kupfer-------------------------------------------------------------------------------------------------------------------------------------------------

theta_cu1, theta_cu2, theta_cu3 = np.array_split(theta_cu,3)
N_Cu1, N_Cu2, N_Cu3 = np.array_split(N_cu,3)

table_header = r'''
\begin{longtable}{S[table-format=2.1] S[table-format=3.0] S[table-format=4.1] S[table-format=2.0] S[table-format=3.1] S[table-format=4.0]}
    \caption{Impulsrate $N$ der Kupferröhre in Abh\"angigkeit des Streuwinkels $\theta$}\\
    \toprule
    {$\theta \; [\si{\degree}]$} & {$\text{Impulse} \; N_{\text{Kupfer}} \; [\si{\per\second}]$} & 
    {$\theta \; [\si{\degree}]$} & {$\text{Impulse} \; N_{\text{Kupfer}} \; [\si{\per\second}]$} & 
    {$\theta \; [\si{\degree}]$} & {$\text{Impulse} \; N_{\text{Kupfer}} \; [\si{\per\second}]$}\\
    \cmidrule(lr{0.5em}){1-2} \cmidrule(lr{0.5em}){3-4} \cmidrule(lr{0.5em}){5-6}
'''
table_footer = r''' 
    \bottomrule
    \label{tab:1}
\end{longtable}
'''

row_template = r'     {0:2.1f} & {1:3.0f} & {2:4.1f} & {3:2.0f} & {4:3.1f} & {5:4.0f}\\'

with open('build/table_Cu.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_cu1, N_Cu1, theta_cu2, N_Cu2, theta_cu3, N_Cu3):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)


#Tabelle Bragg-------------------------------------------------------------------------------------------------------------------------------------------------

theta_bragg1, theta_bragg2 = np.array_split(theta_bragg, 2)
N_bragg1, N_bragg2 = np.array_split(N_bragg, 2)

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    {$\text{Winkel} \; \theta\; [\si{\degree}]$} & {$\text{Impulse $N_{\text{Bragg}}$} \; [\si{\per\second}]$} &
    {$\text{Winkel} \; \theta\; [\si{\degree}]$} & {$\text{Impulse $N_{\text{Bragg}}$} \; [\si{\per\second}]$} \\

    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f} & {2:1.2f} & {3:1.2f}  \\'

with open('build/table_bragg.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_bragg1, N_bragg1, theta_bragg2, N_bragg2):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

#Tabelle Brom-------------------------------------------------------------------------------------------------------------------------------------------------

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    {$\text{Winkel} \; \theta\; [\si{\degree}]$} & {$\text{Impulse $N_{\text{Brom}}$} \; [\si{\per\second}]$} \\

    \cmidrule(lr{0,5em}){1-2}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f}  \\'

with open('build/table_brom.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_brom, N_brom):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

#Tabelle Gallium-------------------------------------------------------------------------------------------------------------------------------------------------

theta_gallium1, theta_gallium2 = np.array_split(theta_gallium, 2)
N_gallium1, N_gallium2 = np.array_split(N_gallium, 2)

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    {$\text{Winkel} \; \theta\; [\si{\degree}]$} & {$\text{Impulse $N_{\text{Gallium}}$} \; [\si{\per\second}]$} &
    {$\text{Winkel} \; \theta\; [\si{\degree}]$} & {$\text{Impulse $N_{\text{Gallium}}$} \; [\si{\per\second}]$} \\

    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f} & {2:1.2f} & {3:1.2f}  \\'

with open('build/table_gallium.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_gallium1, N_gallium1, theta_gallium2, N_gallium2):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

#Tabelle Rubidium-------------------------------------------------------------------------------------------------------------------------------------------------

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    {$\text{Winkel} \; \theta\; [\si{\degree}]$} & {$\text{Impulse $N_{\text{Rubidium}}$} \; [\si{\per\second}]$} \\

    \cmidrule(lr{0,5em}){1-2}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f}  \\'

with open('build/table_rub.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_rub, N_rub):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

#Tabelle Rubidium-------------------------------------------------------------------------------------------------------------------------------------------------

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    {$\text{Winkel} \; \theta\; [\si{\degree}]$} & {$\text{Impulse $N_{\text{Strontium}}$} \; [\si{\per\second}]$} \\

    \cmidrule(lr{0,5em}){1-2}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f}  \\'

with open('build/table_stron.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_stron, N_stron):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

#Tabelle Zink-------------------------------------------------------------------------------------------------------------------------------------------------

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    {$\text{Winkel} \; \theta\; [\si{\degree}]$} & {$\text{Impulse $N_{\text{Zink}}$} \; [\si{\per\second}]$} \\

    \cmidrule(lr{0,5em}){1-2}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f}  \\'

with open('build/table_zink.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_zink, N_zink):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

#Tabelle Zirkonium-------------------------------------------------------------------------------------------------------------------------------------------------

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    {$\text{Winkel} \; \theta\; [\si{\degree}]$} & {$\text{Impulse $N_{\text{Zirkonium}}$} \; [\si{\per\second}]$} \\

    \cmidrule(lr{0,5em}){1-2}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f}  \\'

with open('build/table_zirk.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_zirk, N_zirk):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)