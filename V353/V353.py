
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
import scipy.constants as const #Bsp.: const.physical_constants["proton mass"], output -> value, unit, error

U0 = 0.3 

f, A, t = np.genfromtxt('data.txt', unpack=True) #Variablen definieren f=Frequenz, A=Amplitude, t=Zeit

def f(f, A, x, t):
    return [...] * np.exp(xt)

A0 = A / U0
#params, cov_mat = np.polyfit(f, A0, deg=1, cov=True)
slope, intercept, r_value, p_value, std_err = stats.linregress(f, A0)
plt.plot(f, A0, 'k.', label="Messwerte")
plt.yscale('log')

x_plot = np.linspace(0, 100)



plt.plot(x_plot, y_plot, 'r-', label="Nicht-lineare Regression")

plt.legend(loc="best")
plt.title('Messwerte + Linear Regression')
plt.xlabel('Frequenz in Hertz')
plt.ylabel('Amplitude in Volt')
plt.savefig('build/plot.pdf')