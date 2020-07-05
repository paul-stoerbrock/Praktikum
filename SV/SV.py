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

q =
d =
rho =
g = const.g 
h =
b =
p =
v_f =
v_r =
V =
E =

def a(b, p, eta, v_f, g, rho):
    return np.sqrt( (b/(2*p))**2 + (9*eta*v_f)/(2*g*rho))-b/(2*p)

def m(a, rho):
    return (4/3)*const.pi*(a**3)*rho

def q(m, rho), rho, g, v_f, v_r, E):
    return (m*g*(v_f+V_r))/(E*v_f)




## tex file for U_g_gruen_err 
#
#with open('build/U_g_gruen_err.tex', 'w') as f:
#  f.write(make_SI(U_g_gruen_err,r'\volt', figures=1))