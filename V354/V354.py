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


# Funktionsdefinitionen #########################################################################



# Messdaten ##################################################################################
Aa, ta = np.genfromtxt('4a.txt', unpack=True) # Aa=Amplitudenspannung für 4a), ta=Zeitdifferenz der Nulldurchgänge für 4a)

f, A, t = np.genfromtxt('data.txt', unpack=True) # f=Frequenz, A=Amplitudenspannung, t=Zeitdifferenz der Nulldurhgänge

U0 = 10 # angelegte Spannung in Volt

A0 = A/U0

phi = f * t * 2 * np.pi # Umrechnung von t in phi


#Plot für 4a) #################################################################################

plt.plot(ta*1e+04,np.log(Aa), 'kx', label="Messdaten")
x_plot = np.linspace(0, 2e-04, 10)
plt.plot(x_plot*1e04,intercept+slope*x_plot, 'r-', label="Lineare Regression")
plt.legend(loc="best")
plt.xlabel('Zeit in s')
plt.ylabel('Log(A) in Volt')
plt.tight_layout
plt.title('Plot für Aufg. 4a)')
plt.savefig('build/plota.pdf')
plt.close()


# Erstellung des Plots c in ln-Darstellung

plt.plot(f, A0, 'rx', label="Messdaten")
plt.xscale('log')
plt.xlabel('Frequenz in Hertz')
plt.ylabel('$U_C/U_0$ in Volt')
plt.tight_layout
plt.savefig('build/plotcln.pdf')
plt.close()

# Erstellung des Plots c in linearer Darstellung

plt.plot(f, A0, 'rx', label="Messdaten")
plt.legend(loc="best")
plt.xlabel('Frequenz in Hertz')
plt.ylabel('$U_C/U_0$ in Volt')
plt.tight_layout
plt.savefig('build/plotclin.pdf')
plt.close()


# Plot d) ###########################################################################################

# Erstellung des Plots d) linear

plt.plot(f*1e-04, phi, 'kx', label="Messdaten")
plt.xscale('log')
plt.yscale('log')
plt.axis([0, 5, 0, np.pi])
plt.legend(loc="best")
plt.xlabel(f'$\nu$ / $10^4 Hz$')
plt.ylabel('Phi in Bogenmaß')
plt.tight_layout
plt.savefig('build/plotd.pdf')
plt.close()

# Erstellung Tabelle c), d) ###################################################################################

f1, f2, f3= np.array_split(f*1e-03, 3)
A1, A2, A3 = np.array_split(A, 3)
t1, t2, t3 = np.array_split(t, 3)


table_header = r'''
  \begin{tabular}{c c c c c c c c c}
    \toprule
    {$\nu \:/\: \si{\kilo\hertz}$}& {$U_C \:/\: \si{\volt}$} & {$t \:/\: \si{\micro\second}$} &
    {$\nu \:/\: \si{\kilo\hertz}$}& {$U_C \:/\: \si{\volt}$} & {$t \:/\: \si{\micro\second}$} &
    {$\nu \:/\: \si{\kilo\hertz}$}& {$U_C \:/\: \si{\volt}$} & {$t \:/\: \si{\micro\second}$} \\


    \cmidrule(lr{0,5em}){1-3} \cmidrule(lr{0,5em}){4-6} \cmidrule(lr{0,5em}){7-9}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f} & {2:1.2f} & {0:1.2f} & {1:1.2f} & {2:1.2f} & {0:1.2f} & {1:1.2f} & {2:1.2f}  \\'

# Tabelle für 4a) wird im Tex Format geschrieben ############################################################################################################################################################################

with open('build/table_c.tex', 'w') as h:
    h.write(table_header)
    for row in zip(f1, A1, t1, f2, A2, t2, f3, A3, t3):
        h.write(row_template.format(*row))
        h.write('\n')
    h.write(table_footer)



# Kontrollprints

