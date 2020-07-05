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

def make_SI(num, unit, exp='', figures=None):
    ''' Format an uncertainties ufloat as a \SI quantity '''
    if np.any(stds([num])):
        if figures is None:
            figures = ''
        x = '{0:.{1:}uf}'.format(num, figures).replace('/', '')
    else:
        x = '{0:.{1:}f}'.format(num, figures)

    return r'\SI{{{}{}}}{{{}}}'.format(x, exp, unit)

t_auf_1, x_auf_1, t_ab_1, x_ab_1 = np.genfromtxt('tropfen_1.txt', unpack=True)
t_auf_2, x_auf_2, t_ab_2, x_ab_2 = np.genfromtxt('tropfen_2.txt', unpack=True)
t_auf_3, x_auf_3, t_ab_3, x_ab_3 = np.genfromtxt('tropfen_3.txt', unpack=True)
t_auf_4, x_auf_4, t_ab_4, x_ab_4 = np.genfromtxt('tropfen_4.txt', unpack=True)
t_auf_5, x_auf_5, t_ab_5, x_ab_5 = np.genfromtxt('tropfen_5.txt', unpack=True)
t_auf_6, x_auf_6, t_ab_6, x_ab_6 = np.genfromtxt('tropfen_6.txt', unpack=True)
t_auf_7, x_auf_7, t_ab_7, x_ab_7 = np.genfromtxt('tropfen_7.txt', unpack=True)

d = 7.6e-3 #in m
rho = 886 #kg/m^3
g = const.g 
eta = 1.843*10**(-5) #Nsm^(-2)
b = 8.2*10**(-3) #in Pa*m
p = 10**5 #in pascal
E = 308/d #in V/m

#Geschwindigkeit Aufwärts============================================================================================================================================================================================================================

v_auf_1_n = (x_auf_1*1e-3)/t_auf_1
v_auf_2_n = (x_auf_2*1e-3)/t_auf_2
v_auf_3_n = (x_auf_3*1e-3)/t_auf_3
v_auf_4_n = (x_auf_4*1e-3)/t_auf_4
v_auf_5_n = (x_auf_5*1e-3)/t_auf_5
v_auf_6_n = (x_auf_6*1e-3)/t_auf_6
v_auf_7_n = (x_auf_7*1e-3)/t_auf_7

v_auf_1_err = sem(v_auf_1_n)
v_auf_2_err = sem(v_auf_2_n)
v_auf_3_err = sem(v_auf_3_n)
v_auf_4_err = sem(v_auf_4_n)
v_auf_5_err = sem(v_auf_5_n)
v_auf_6_err = sem(v_auf_6_n)
v_auf_7_err = sem(v_auf_7_n)

v_r_n = np.array([np.mean(v_auf_1_n), np.mean(v_auf_2_n), np.mean(v_auf_3_n), np.mean(v_auf_4_n), np.mean(v_auf_5_n), np.mean(v_auf_6_n), np.mean(v_auf_7_n)])
v_r_err = np.array([v_auf_1_err, v_auf_2_err, v_auf_3_err, v_auf_4_err, v_auf_5_err, v_auf_6_err, v_auf_7_err])
v_r = unp.uarray([v_r_n, v_r_err])

#Geschwindigkeit Abwärts============================================================================================================================================================================================================================

v_ab_1_n = (x_ab_1*1e-3)/t_ab_1
v_ab_2_n = (x_ab_2*1e-3)/t_ab_2
v_ab_3_n = (x_ab_3*1e-3)/t_ab_3
v_ab_4_n = (x_ab_4*1e-3)/t_ab_4
v_ab_5_n = (x_ab_5*1e-3)/t_ab_5
v_ab_6_n = (x_ab_6*1e-3)/t_ab_6
v_ab_7_n = (x_ab_7*1e-3)/t_ab_7

v_ab_1_err = sem(v_ab_1_n)
v_ab_2_err = sem(v_ab_2_n)
v_ab_3_err = sem(v_ab_3_n)
v_ab_4_err = sem(v_ab_4_n)
v_ab_5_err = sem(v_ab_5_n)
v_ab_6_err = sem(v_ab_6_n)
v_ab_7_err = sem(v_ab_7_n)

v_f_n = np.array([np.mean(v_ab_1_n), np.mean(v_ab_2_n), np.mean(v_ab_3_n), np.mean(v_ab_4_n), np.mean(v_ab_5_n), np.mean(v_ab_6_n), np.mean(v_ab_7_n)])
v_f_err = np.array([v_ab_1_err, v_ab_2_err, v_ab_3_err, v_ab_4_err, v_ab_5_err, v_ab_6_err, v_ab_7_err])
v_f = unp.uarray([v_f_n, v_f_err])

#Funktionsdefinitionen============================================================================================================================================================================================================================

def a(b, p, eta, v_f, g, rho):
    return unp.sqrt( (b/(2*p))**2 + (9*eta*v_f)/(2*g*rho))-b/(2*p)

def m(a, rho):
    return (4/3)*const.pi*(a**3)*rho

def q(m, rho, g, v_f, v_r, E):
    return (m*g*(v_f+v_r))/(E*v_f)

print(a(b, p, eta, v_f, g, rho))
print(q(m(a(b, p, eta, v_f, g, rho), rho), rho, g, v_f, v_r, E)/const.e)

#Plot============================================================================================================================================================================================================================

plt.plot(noms(a(b, p, eta, v_f, g, rho))*1e+3, noms(q(m(a(b, p, eta, v_f, g, rho), rho), rho, g, v_f, v_r, E))/const.e, 'kx', label='Messwerte')

plt.xlabel(r'Radius des Tropfens $\;[mm]$')
plt.ylabel(r'Vielfaches der Elementarladung')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('plot.pdf')
plt.close()