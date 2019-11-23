
import pylab
import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
from uncertainties.unumpy import (
    nominal_values as noms,
    std_devs as stds
)
from scipy import stats 
from scipy.stats import sem #standard error of mean = sem(x)
from scipy.optimize import curve_fit #function curve_fit 
import scipy.constants as const #Bsp.: const.physical_constants["proton mass"], output -> value, unit, error

# A(Omega)/U0 ############################################################################################################################################################################
def UC(U, tc, RC):
    return np.log(U0 (1 - U**((-1/RC)* tc)))

def g(f, RC):
    return  1/(np.sqrt(4 * (np.pi)**2 * f**2 * RC**2+1))

def h(f, RC):
    return np.arctan(-2 * np.pi * f * RC)

def d(A0 ,f , RC):
    return np.arcsin(-2 * np.pi * f * RC * A0)


U0 = 3 #kann es sein, dass U0 eigentlich 3V war, und nicht 0,3V?

U, tc = np.genfromtxt('4a.txt', unpack=True) #Variablen definieren U=Q/C, tc=Zeit aus 4a)

f, A, t = np.genfromtxt('data.txt', unpack=True) #Variablen definieren f=Frequenz, A=Amplitude, t=Zeit
t *= 1e-03
A0 = A / U0
phi = f * t * 2 * np.pi

# Definition von RC (popt) ############################################################################################################################################################################

slope, intercept, r_value, p_value, std_err = stats.linregress(tc, U)

popt, pcov = curve_fit(
    g,
    f,
    A0,
    sigma=None,
    absolute_sigma=True,
    p0=[1e-03]
    )

phiRC, phicov = curve_fit(
    h,
    f,
    phi,
    sigma=None,
    absolute_sigma=True,
    p0=[1e-03]
    )


# Plot für 4b) [A(w)/U0] ############################################################################################################################################################################

plt.plot(tc, U, 'kx', label="4a) RC aus Ausgleichsgerade")
plt.xscale('log')
plt.plot(tc, intercept + slope*tc, 'r-', label="Lineare Regression")
plt.legend(loc="best")
plt.title('4a)')
plt.xlabel('Zeit in ms')
plt.ylabel('$U_c$ in Volt')
plt.tight_layout
plt.savefig('build/4a.pdf')
plt.close()

plt.plot(f, A0, 'kx', label="Frequenz und Amplitude")
plt.xscale('log')
x_plot = np.linspace(1, 100000, 100000)
plt.plot(x_plot, g(x_plot,*popt), 'r-', label="Nicht-lineare Regression")
plt.legend(loc="best")
plt.title('4b)')
plt.xlabel('Frequenz in Hertz')
plt.ylabel('A/$U_0$ in Volt')
plt.tight_layout
plt.savefig('build/plotA.pdf')
plt.close()

plt.plot(f, phi, 'kx', label="Frequenz und Phase")
plt.xscale('log')
x_plot = np.linspace(1, 1000, 1000)
plt.plot(x_plot, h(x_plot,*phiRC), 'r-', label="Nicht-lineare Regression")
plt.legend(loc="best")
plt.title('4c)')
plt.xlabel('Frequenz in Hertz')
plt.ylabel('Phase')
plt.grid(True)
plt.tight_layout
plt.savefig('build/plot.pdf')
plt.close()

plt.polar(phi, A0, 'kx', label="Amplitude und Phase")
x_plot = np.linspace(0, 1, 24)
plt.polar(d(x_plot,f , phiRC), x_plot, 'r-', label="Polarplot")
plt.legend(loc="best")
#plt.title('4d)')
plt.grid(True)
plt.tight_layout
plt.savefig('build/plotpolar.pdf')



#SI Einheiten

ohmF=r' }{\Ohm\Farad}$'

with open('build/mean_RC.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{popt[0]:.2f}')
    RC.write(ohmF)

with open('build/mean_phiRC.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{phiRC[0]:.2f}')
    RC.write(ohmF)

# Tabelle für 4a) wird erstellt ############################################################################################################################################################################

U1, U2, U3= np.array_split(U, 3)
tc1, tc2, tc3 = np.array_split(tc, 3)

table_header = r'''
  \begin{tabular}{c c c c c c}
    \toprule
    {$U_C \:/\: \si{\volt}$} & {$t \:/\: \si{\milli\second}$} & 
    {$U_C \:/\: \si{\volt}$} & {$t \:/\: \si{\milli\second}$} &  
    {$U_C \:/\: \si{\volt}$} & {$t \:/\: \si{\milli\second}$}\\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4} \cmidrule(lr{0,5em}){5-6}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f} & {2:1.2f} & {3:1.2f} & {4:1.2f} & {5:1.2f}\\'

# Tabelle für 4a) wird im Tex Format geschrieben ############################################################################################################################################################################

with open('build/table_4a.tex', 'w') as h:
    h.write(table_header)
    for row in zip(U1, tc1, U2, tc2, U3, tc3):
        h.write(row_template.format(*row))
        h.write('\n')
    h.write(table_footer)

# Tabelle für 4b-c) wird erstellt ############################################################################################################################################################################

f1, f2, f3= np.array_split(f, 3)
A1, A2, A3 = np.array_split(A, 3)
t1, t2, t3 = np.array_split(1000*t, 3)

table_header = r'''
  \begin{tabular}{c c c c c c c c c}
    \toprule
    {$f \:/\: \si{\hertz}$} & {$A(\omega) \:/\: \si{\milli\volt}$} & {$\Delta T \:/\: \si{\micro\second}$} & 
    {$f \:/\: \si{\hertz}$} & {$A(\omega) \:/\: \si{\milli\volt}$} & {$\Delta T \:/\: \si{\micro\second}$} & 
    {$f \:/\: \si{\hertz}$} & {$A(\omega) \:/\: \si{\milli\volt}$} & {$\Delta T \:/\: \si{\micro\second}$}\\
    \cmidrule(lr{0,5em}){1-3} \cmidrule(lr{0,5em}){4-6} \cmidrule(lr{0,5em}){7-9}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.0f} & {1:1.2f} & {2:1.2f} & {3:1.0f} & {4:1.2f} & {5:1.2f} & {6:1.0f} & {7:1.2f} & {8:1.2f} \\'

# Tabelle für 4b-c) wird im Tex Format geschrieben ############################################################################################################################################################################

with open('build/table_4b.tex', 'w') as i:
    i.write(table_header)
    for row in zip(f1, A1*1e+03, t1, f2, A2*1e+03, t2, f3, A3*1e+03, t3):
        i.write(row_template.format(*row))
        i.write('\n')
    i.write(table_footer)


# Testprints #########################################################################

print(-(np.sin(phi))/(2 * np.pi * f * phiRC[0]))
print(phi)
print(A0)