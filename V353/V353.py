
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



U0 = 3 #kann es sein, dass U0 eigentlich 3V war, und nicht 0,3V?

f, A, t = np.genfromtxt('data.txt', unpack=True) #Variablen definieren f=Frequenz, A=Amplitude, t=Zeit

#A(Omega)/U0

def g(f, RC):
    return  1/(np.sqrt(4 * (np.pi)**2 * f**2 * RC**2+1))

A0 = A / U0

#Definition von RC (popt)

popt, pcov = curve_fit(
    g,
    f,
    A0,
    sigma=None,
    absolute_sigma=True,
    p0=[1e-03]
    )
print(popt) #Überprüfung von Größe RC

plt.plot(f, A0, 'kx', label="Messwerte")
plt.xscale('log')

x_plot = np.linspace(1, 100000, 100000)

plt.plot(x_plot, g(x_plot,*popt), 'r-', label="Nicht-lineare Regression")

plt.legend(loc="best")
plt.title('Messwerte + non-linear regression')
plt.xlabel('Frequenz in Hertz')
plt.ylabel('A/$U_0$ in Volt')
plt.grid(True)
plt.tight_layout
plt.savefig('build/plot.pdf')

#SI Einheiten

mean_RC=r' }{\Ohm\Farad}$'

with open('build/mean_RC.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{popt[0]:.2f}')
    RC.write(mean_RC)

#Tabelle wird erstellt

table_header = r'''
  \begin{tabular}{c c c}
    \toprule
    {$f \:/\: \si{\hertz}$} & {$A(\omega) \:/\: \si{\milli\volt}$} & {$\Delta T \:/\: \si{\milli\second}$}\\
    \midrule
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.0f} & {1:1.2f} & {2:1.2f}  \\'

#Tabelle wird im Tex Format geschrieben

with open('build/table_4b.tex', 'w') as g:
    g.write(table_header)
    for row in zip(f, A*1e+03, t*1e+03):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)