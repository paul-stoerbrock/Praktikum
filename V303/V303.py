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

def U(phi, U0, B):
    return 2/np.pi* U0 * np.cos(phi+B)

def U2(phi, A, B, C, D):
    return A * np.cos(B*phi+C)+D

def y(r, A, B):
    return A*1/r**2+B

# Definition der Messdaten ##############################################################################################

# Messwerte ohne Störung
U_oS, phi_oS = np.genfromtxt('aufg2.txt', unpack=True)# U_oS = Spannung ohne Störung in Volt, phi_oS = Phasenverschiebung ohne Störung in Gradmaß

phi_oS_rad = phi_oS *2*np.pi/360

# Messwerte mit Störung
U_mS, phi_mS = np.genfromtxt('aufg3.txt', unpack=True)# U_oS = Spannung mit Störung in Volt, phi_oS = Phasenverschiebung mit Störung in Gradmaß

phi_mS_rad = phi_mS *2*np.pi/360

# Messwerte der Photodetektorschaltung
r, I = np.genfromtxt('aufg4.txt', unpack=True)# r = Abstand zur Quelle in cm, I = Intensität des Signals in Volt

r*=1e-02

# Plots ################################################################################################################

# Plot für den Fall ohne Störung

par, covm = curve_fit(
    U,
    phi_oS_rad,
    U_oS,
    sigma=None,
    absolute_sigma=True,
    p0=[10, 1]
    )
err = np.sqrt(np.diag(covm))


plt.plot(phi_oS_rad, U_oS,'kx', label='Messwerte')
x_plot = np.linspace(0, 2*np.pi, 1000)
plt.plot(x_plot,U(x_plot, *par) , 'r-', label="Nicht-Lineare Regression")
plt.legend(loc="best")
plt.xlabel(r'Phase $\phi\:/\:rad$')
plt.ylabel(r'Spannung $U\:/\:V$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotphi_oS.pdf')
plt.close()

parphi_oS=unp.uarray(par, err)

print(par)

# Plot für den Fall mit Störung

par, covm = curve_fit(
    U,
    phi_mS_rad,
    U_mS,
    sigma=None,
    absolute_sigma=True,
    p0=[10, 1]
    )
err = np.sqrt(np.diag(covm))


plt.plot(phi_mS_rad, U_mS,'kx', label='Messwerte')
x_plot = np.linspace(0, 2*np.pi, 1000)
plt.plot(x_plot,U(x_plot, *par) , 'r-', label="Nicht-Lineare Regression")
plt.legend(loc="best")
plt.xlabel(r'Phase $\phi\:/\:rad$')
plt.ylabel(r'Spannung $\:/\:V$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotphi_mS.pdf')
plt.close()

parphi_mS=unp.uarray(par, err)

print(par)

# Plot für den Photodetektorschaltung
par, covm = curve_fit(
    y,
    r[1:],
    I[1:],
    sigma=None,
    absolute_sigma=True,
    p0=[100,30]
    )
err = np.sqrt(np.diag(covm))

plt.plot(r[1:], I[1:],'kx', label='Messwerte')
x_plot = np.linspace(0.01, 0.5, 100)
plt.plot(x_plot, y(x_plot, *par))
plt.legend(loc="best")
plt.xlabel(r'Radius$\:r\:/\:m$')
plt.ylabel(r'Intensität$\:U\:/\:V$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotI.pdf')
plt.close()

parI=unp.uarray(par, err)

print(par)

# .tex ####################################################################################################################

# tex file for A_os

with open('build/A_os.tex', 'w') as f:
  f.write(make_SI(parphi_oS[0],r'', figures=2))

# tex file for B_os

with open('build/B_os.tex', 'w') as f:
  f.write(make_SI(parphi_oS[1],r'', figures=2))

# tex file for A_ms

with open('build/A_ms.tex', 'w') as f:
  f.write(make_SI(parphi_mS[0],r'', figures=2))

# tex file for B_ms

with open('build/B_ms.tex', 'w') as f:
  f.write(make_SI(parphi_mS[1],r'', figures=2))

# tex file for I_A

with open('build/I_A.tex', 'w') as f:
  f.write(make_SI(parI[0]*1e03,r'',exp='e-03' ,figures=2))

# tex file for I_B

with open('build/I_B.tex', 'w') as f:
  f.write(make_SI(parI[1],r'', figures=2))



# Tabellen ################################################################################################################

# Ohne Störung =================================================================================================================

U_oS1, U_oS2 = np.array_split(U_oS,2)
phi_oS1, phi_oS2 = np.array_split(phi_oS,2)

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    \multicolumn{1}{c}{$\text{Spannung} \; U\:/\:\si{\volt}$} & \multicolumn{1}{c}{$\text{Phase} \; \varphi\:/\:\si{\degree}$} & 
    \multicolumn{1}{c}{$\text{Spannung} \; U\:/\:\si{\volt}$} & \multicolumn{1}{c}{$\text{Phase} \; \varphi\:/\:\si{\degree}$}\\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4}

'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.2f} & {2:1.2f} & {3:1.2f}\\'

with open('build/table_oS.tex', 'w') as g:
    g.write(table_header)
    for row in zip(U_oS1, phi_oS1, U_oS2, phi_oS2):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

# Mit Störung =================================================================================================================

U_mS1, U_mS2 = np.array_split(U_mS,2)
phi_mS1, phi_mS2 = np.array_split(phi_mS,2)

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    \multicolumn{1}{c}{$\text{Spannung} \; U\:/\:\si{\volt}$} & \multicolumn{1}{c}{$\text{Phase} \; \varphi\:/\:\si{\degree}$} & 
    \multicolumn{1}{c}{$\text{Spannung} \; U\:/\:\si{\volt}$} & \multicolumn{1}{c}{$\text{Phase} \; \varphi\:/\:\si{\degree}$}\\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.2f} & {2:1.2f} & {3:1.2f}\\'

with open('build/table_mS.tex', 'w') as g:
    g.write(table_header)
    for row in zip(U_mS1, phi_mS1, U_mS2, phi_mS2):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

# Intensität =================================================================================================================

r1, r2 = np.array_split(r,2)
I1, I2 = np.array_split(I,2)

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    \multicolumn{1}{c}{$\text{Radius} \; r\:/\si{\centi\meter}$} & \multicolumn{1}{c}{$\text{Intensität} \; I\:/\:\si{\volt} $} & 
    \multicolumn{1}{c}{$\text{Radius} \; r\:/\si{\centi\meter}$} & \multicolumn{1}{c}{$\text{Intensität} \; I\:/\:\si{\volt}$}\\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.2f} & {2:1.2f} & {3:1.2f}\\'

with open('build/table_I.tex', 'w') as g:
    g.write(table_header)
    for row in zip(r1*1e02, I1, r2*1e02, I2):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)