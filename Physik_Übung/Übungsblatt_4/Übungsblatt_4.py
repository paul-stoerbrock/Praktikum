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


def f_real(E, V_2):
     return np.cos(2*(-5*np.sqrt(2*E)+np.arctan(np.sqrt(E/(V_2-E)))+np.pi-np.arctan(np.sqrt((E)/(10-E)))))-1

def f_imag(E, V_2):
     return np.sin(2*(-5*np.sqrt(2*E)+np.arctan(np.sqrt(E/(V_2-E)))+np.pi-np.arctan(np.sqrt((E)/(10-E)))))


c=9

x_plot = np.linspace(0, c-0.1, 1000)
plt.plot(x_plot, f_real(x_plot,c),'r-', label="$E$")
plt.plot(x_plot, f_imag(x_plot,c) ,'b.', label="$E$")
plt.legend(loc="best")
#plt.xlabel(r'Energie ')
#plt.ylabel(r'$k_1/k_2$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot.pdf')
plt.close()
