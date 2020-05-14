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

#x_plot = np.linspace(100, 1000, 100000)
#plt.plot(x_plot, np.sqrt(2*x_plot)/(1.05 *1e-34) ,'r-', label="$k_1$")
#plt.plot(x_plot, (np.sqrt(2*(x_plot-100)))/(1.05 *1e-34) ,'b-', label="$k_2$")
#plt.legend(loc="best")
#plt.xlabel(r'Energie ')
#plt.ylabel(r'$k_1/k_2$')
#plt.grid()
#plt.tight_layout
#plt.savefig('build/plot.pdf')
#plt.close()