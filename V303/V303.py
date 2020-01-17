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

# Funktionsdefinitionen ########################################################################################################################################################

def make_SI(num, unit, exp='', figures=None):
    ''' Format an uncertainties ufloat as a \SI quantity '''
    if np.any(stds([num])):
        if figures is None:
            figures = ''
        x = '{0:.{1:}uf}'.format(num, figures).replace('/', '')
    else:
        x = '{0:.{1:}f}'.format(num , figures)

    return r'\SI{{{}{}}}{{{}}}'.format(x, exp, unit)

def U(phi, U0, B):
    return 2/np.pi* U0 * np.cos(phi+B)

def y(r, A):
    return A*1/r**2

# Definition der Messdaten ##############################################################################################

# Messwerte ohne Störung
U_oS, phi_oS = np.genfromtxt('aufg2.txt', unpack=True)# U_oS = Spannung ohne Störung in Volt, phi_oS = Phasenverschiebung ohne Störung in Gradmaß

phi_oS_rad = phi_oS *2*np.pi/360

# Messwerte mit Störung
U_mS, phi_mS = np.genfromtxt('aufg3.txt', unpack=True)# U_oS = Spannung mit Störung in Volt, phi_oS = Phasenverschiebung mit Störung in Gradmaß

phi_mS_rad = phi_mS *2*np.pi/360

# Messwerte der Photodetektorschaltung
r, I = np.genfromtxt('aufg4.txt', unpack=True)# r = Abstand zur Quelle in cm, I = Intensität des Signals in Volt

r*=1e-02

# Plots ################################################################################################################

# Plot für den Fall ohne Störung

par, covm = curve_fit(
    U,
    phi_oS_rad,
    U_oS,
    sigma=None,
    absolute_sigma=True,
    p0=[10, 1]
    )
err = np.sqrt(np.diag(covm))


plt.plot(phi_oS_rad, U_oS,'kx', label='Messwerte')
x_plot = np.linspace(0, 2*np.pi, 1000)
plt.plot(x_plot,U(x_plot, *par) , 'r-', label="Nicht-Lineare Regression")
plt.legend(loc="best")
plt.xlabel(r'Phase $\phi\:/\:rad$')
plt.ylabel(r'Spannung $U\:/\:V$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotphi_oS.pdf')
plt.close()

parphi_oS=unp.uarray(par, err)

print(par)

# Plot für den Fall mit Störung

par, covm = curve_fit(
    U,
    phi_mS_rad,
    U_mS,
    sigma=None,
    absolute_sigma=True,
    p0=[10, 1]
    )
err = np.sqrt(np.diag(covm))


plt.plot(phi_mS_rad, U_mS,'kx', label='Messwerte')
x_plot = np.linspace(0, 2*np.pi, 1000)
plt.plot(x_plot,U(x_plot, *par) , 'r-', label="Nicht-Lineare Regression")
plt.legend(loc="best")
plt.xlabel(r'Phase $\phi\:/\:rad$')
plt.ylabel(r'Spannung $\:/\:V$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotphi_mS.pdf')
plt.close()

parphi_mS=unp.uarray(par, err)

print(par)

# Plot für den Photodetektorschaltung
par, covm = curve_fit(
    y,
    r[1:],
    I[1:],
    sigma=None,
    absolute_sigma=True,
    p0=[1]
    )
err = np.sqrt(np.diag(covm))

plt.plot(r[1:], I[1:],'kx', label='Messwerte')
x_plot = np.linspace(0.01, 0.5, 100)
plt.plot(x_plot, y(x_plot, *par))
plt.legend(loc="best")
plt.xlabel(r'Radius$\:r\:/\:m$')
plt.ylabel(r'Intensität\:$U\:/\:V$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotI.pdf')
plt.close()

parI=unp.uarray(par, err)

print(par)