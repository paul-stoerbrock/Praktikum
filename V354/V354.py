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

#def U_C(f, L, R, C):
#    return 1/(np.sqrt(R**2 + (2*np.pi*f*L-1/(2*np.pi*f*C))**2))


# Messdaten ##################################################################################

f, A, t = np.genfromtxt('data.txt', unpack=True) # f=Frequenz, A=Amplitudenspannung, t=Zeitdifferenz der Nulldurhgänge

U0 = 10 # angelegte Spannung in Volt

A0 = A/U0

t *= 1e+06

phi = f * t * 2 * np.pi # Umrechnung von t in phi




#Plot für c) #################################################################################

# Erstellung des Plots c in ln-Darstellung

plt.plot(f, A0, 'rx', label="Messdaten")
plt.xscale('log')
plt.legend(loc="best")
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

