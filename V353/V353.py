
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

U0 = 0.3 

f, A, t = np.genfromtxt('data.txt', unpack=True) #Variablen definieren f=Frequenz, A=Amplitude, t=Zeit

def g(f, RC):
    return  1/(np.sqrt(4 * (np.pi)**2 * f**2 * RC**2+1))

A0 = A / U0


popt, pcov = curve_fit(
    g,
    f,
    A0,
    sigma=None,
    absolute_sigma=True,
    p0=[1e-03]
    )
print(popt) #Überprüfung von Größe RC

plt.plot(f, A0, 'k.', label="Messwerte")
plt.yscale('log')

x_plot = np.linspace(0, 100000, 100000)



plt.plot(x_plot, g(x_plot,*popt), 'r-', label="Nicht-lineare Regression")

plt.legend(loc="best")
plt.title('Messwerte + Linear Regression')
plt.xlabel('Frequenz in Hertz')
plt.ylabel('A/U0 in Volt')
plt.tight_layout
plt.savefig('build/plot.pdf')