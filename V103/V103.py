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

# Funktionsdefinitionen

def make_SI(num, unit, exp='', figures=None):
    ''' Format an uncertainties ufloat as a \SI quantity '''
    if np.any(stds([num])):
        if figures is None:
            figures = ''
        x = '{0:.{1:}uf}'.format(num, figures).replace('/', '')
    else:
        x = '{0:.{1:}f}'.format(num, figures)

    return r'\SI{{{}{}}}{{{}}}'.format(x, exp, unit)

# Messwerte ####################################################################################################################

xCu_einohne, DCu_einohne,xCu_einmit, DCu_einmit = np.genfromtxt('dataCuein.txt', unpack=True) # Messwette für Kupfer einseitig belastet 
xCu_dopohne, DCu_dopohne,xCu_dopmit, DCu_dopmit = np.genfromtxt('dataCudop.txt', unpack=True) # Messwerte für Kupfer doppelseitig belastet
xAl_einohne, DAL_einohne,xAl_einmit, DAl_einmit = np.genfromtxt('dataAlein.txt', unpack=True) # Messwerte für Aluminium einseitig belastet
xAl_dopohne, DAl_dopohne,xAl_dopmit, DAl_dopmit = np.genfromtxt('dataAldop.txt', unpack=True) # Messwerte für Aluminium doppelseitig belastet

# Maße der Stäbe
l_CU = 0.600
l_Al = 0.592
r_Al =np.array([10.00, 10.00, 10.10, 10.00, 10.10, 10.10, 10.10, 10.00, 10.25, 10.00])/2*1e-03
d_CU =np.array([10.02, 10.05, 10.06, 10.03, 10.04, 10.09, 10.04, 10.05, 10.08, 10.05])*1e-03
b_CU =np.array([10.06, 10.02, 10.05, 10.01, 10.03, 10.03, 10.04, 10.02, 10.04, 10.03])*1e-03

# Massen der Körper
m_aufhaeng = 19 *1e-03
m_schraube = 22.1*1e-03
m_Cu_stange = 528.8*1e-03
m_Al_stange = 132.6*1e-03

# Gewichte Cuein
m_Cuein1 = 502.5*1e-03
m_Cuein2 = 503.3*1e-03

# Gewichte Cudop
m_Cudop1 = 1170.5*1e-03
m_Cudop2 = 1159.7*1e-03
m_Cudop3 = 1162.3*1e-03

# Gewichte Alein
m_Alein1 = 500.1*1e-03

# Gewichte Aldop
m_Aldop1 = 500.1*1e-03
m_Aldop2 = 499.8*1e-03
m_Aldop3 = 226.7*1e-03

# Berechnungen ###########################################################################################################

# Berechnung des Flächenträgheitsmoment

# Für Alein




# Erstellung der Plots ###############################################################################################################


