
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

f, A, a = np.genfromtxt('data.txt', unpack=True) #Variablen definieren f=Frequenz, A=Amplitude, a=Zeit

params, cov_mat = np.polyfit(a, A, deg=1, cov=True)
a.set_xscale('linear')
A.set_yscale('log')
plt.plot(a, A, 'k.', label="Messwerte")

x_plot = np.linspace(0, 10)
plt.plot(x_plot, params[0] * x_plot + params[1], 'r-', label="Lineare Regression")
plt.savefig('build/plot.pdf')
plt.legend(loc="best")
plt.title('Messwerte + Linear Regression')
plt.xlabel('Zeit in Sekunden')
plt.ylabel('Amplitude in Volt')
plt.savefig('build/4a.pdf')