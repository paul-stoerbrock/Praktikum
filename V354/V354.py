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


# Messdaten ##################################################################################

f, A, t = np.genfromtxt('data.txt', unpack=True) # f=Frequenz, A=Amplitudenspannung, t=Zeitdifferenz der Nulldurhgänge

t *= 1e+06


U0 = 10 # angelegte Spannung

phi = f * t * 2 * np.pi # Umrechnung von t in phi




#Plot für 4b) #################################################################################

plt.plot(f, A/U0, 'kx', label="Messdaten")
plt.yscale('log')
#plt.plot(tc, np.exp(intercept + slope*tc), 'r-', label="Lineare Regression")
plt.legend(loc="best")
plt.xlabel('Frequenz in Hertz')
plt.ylabel('$U_c/U_0$ in Volt')
plt.tight_layout
plt.savefig('build/plotc.pdf')
plt.close()