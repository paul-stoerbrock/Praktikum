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

# Funktion zur Berechnung des spezifischen Widerstands

def rho(R, l, A):
  return R*A/l

# Funktion zur Berechnung der Ladungsträgerzahl pro Volumen

def n(m, d, c):
  return  -c/(const.e*d*m)

# Funktion zur Berechnung der Hall-Konstante

def AH(n):
  return -1/(n*const.e)

# Funktion zur Berechnung der Ladungsträgerzahl pro Atom

def z(n, varrho, atom_mass):
  return n/(varrho/(atom_mass*const.u))

# Funktion zur Berechnung der mittleren Flugzeit

def tau(n, rho):
  return 2*const.m_e/(const.e**2*n*rho)

# Funktion zur Berechnung der mittleren Driftgeschwindigkeit

def v_d(n):
  return -1e06/(n+const.e)

# Funktion zur Berechnung der Totalgeschwindigkeit

def v(n):
  return unp.sqrt((2*const.h**2/(2*const.m_e)*((3/(8 *np.pi)*n)**2)**(1/3))/const.m_e)

# Funktion zur Berechnung der mittleren freien Weglänge

def l(tau, v):
  return tau*v

# Funktion zur Berechnung der Beweglichkeit

def mu(v_d, n, tau):
  return 1/2*v_d*(const.e**2*n*tau)/(1*10**6*const.m_e)

# Funktion zur Berechnung des B-Feldes

def B1(auf, ab, I):
  return (auf+ab)/2 *I

# Konstanten ########################################################################################################################################################

# Hall- Konstanten [Einheit: m^3 * C^(-1)]:

Hconst_cu_lit = -5.2e-11
Hconst_zn_lit = +6.4e-11
Hconst_ag_lit = -8.9e-11

# Spezifischer Widerstand [Einheit: Ohm * m]:

SpWi_cu_lit = 0.018*1e-06
SpWi_zn_lit = 0.06*1e-06
SpWi_ag_lit = 0.016*1e-06

# Variablen ########################################################################################################################################################

Cu_IB, Cu_UHB, Cu_IQ, Cu_UHQ = np.genfromtxt('cu.txt', unpack=True)

Zn_IB, Zn_UHB, Zn_IQ, Zn_UHQ = np.genfromtxt('zn.txt', unpack=True)

Ag_IB, Ag_UHB, Ag_IQ, Ag_UHQ = np.genfromtxt('ag.txt', unpack=True)

hy_Iauf, hy_Bauf, hy_Iab, hy_Bab = np.genfromtxt('hysterese.txt', unpack=True)

# Werte für Silber %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d_Ag_Folie = 0.205e-03 # Dicke der Silber Folie in m

d_Ag_Draht = 0.205 *1e-03 # Durchmesser des Silberdrahts in m

R_Ag = 0.703 # Widerstand von Aluminiumspule in Ohm

l_Ag = 173 *1e-02 # Länge der Drahtspule in m

atom_mass_Ag = 107.8682 # Atommasse von Silber

varrho_Ag = 10500 # Dichte von Silber in kg/m^3

# Werte für Kupfer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d_Cu_Folie = 18e-06 # Dicke der Kupfer Folie in m

d_Cu_Draht = 0.1e-03 # Durchhmesser des Kupferdrahts in m

R_Cu = 2.903 # Widerstand der Kupferspule in Ohm

l_Cu = 137e-02 # Länge des Kupferdrahtes in m

varrho_Cu = 8960 # Literaturwert für Dichte von Kupfer in kg/m^3

atom_mass_Cu = 63.546 # Atommase von Kupfer in u

# Werte für Zink %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d_Zn_Folie = 0.027e-03 # Dicke der Zink Folie in m

atom_mass_Zn = 65.39 # Atommasse von Zink in u

varrho_Zn = 7130     # Dichte von Zink in kg/m^3 

# Plots ########################################################################################################################################################

#Plot der Hysteresekurve von der Spule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par1, covm1 = np.polyfit(hy_Iauf, hy_Bauf, deg=1, cov=True)
err1 = np.sqrt(np.diag(covm1))

par2, covm2 = np.polyfit(hy_Iab[1:9], hy_Bab[1:9], deg=1, cov=True)
err2 = np.sqrt(np.diag(covm2))


plt.plot(hy_Iauf, hy_Bauf,'kx', label='Messwerte auf')
plt.plot(hy_Iab, hy_Bab,'bx', label='Messwerte ab')
x_plot = np.linspace(0, 5, 1000)
plt.plot(x_plot, x_plot*par1[0]+par1[1], 'r-', label="Lineare Regression Hysteresekurve auf")
plt.plot(x_plot, x_plot*par2[0]+par2[1], 'g-', label="Lineare Regression Hysteresekurve ab")

plt.legend(loc="best")
plt.xlabel(r'Stromstärke $I\:/\:V$')
plt.ylabel(r'magnetische Feldstärke $B\:/\:T$')
plt.grid()
plt.tight_layout
plt.savefig('build/plothy.pdf')
plt.close()

parhyauf=unp.uarray(par1, err1)
parhyab=unp.uarray(par2, err2)

#Plot von Kupfer mit konstantem Querstrom %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B= (parhyauf[0].n + parhyab[0].n)/2*Cu_IB

par, covm = np.polyfit(B, Cu_UHB, deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(B, Cu_UHB,'kx', label='Messwerte')
x_plot = np.linspace(0, 1.3, 1000)
plt.plot(x_plot, x_plot*par[0]+par[1], 'r-', label="Lineare Regression")
plt.yticks([0, 5*1e-06, 10*1e-06, 15*1e-06, 20*1e-06],
           [0, 5, 10, 15, 20])
plt.legend(loc="best")
plt.xlabel(r'magnetische Feldstärke $B\:/\:T$')
plt.ylabel(r'Hallspannung $U_H\:/\:\mu V$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotCu_IB.pdf')
plt.close()

parCu_IB=unp.uarray(par, err)

#Plot von Kupfer mit konstantem B-Feld %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par, covm = np.polyfit(Cu_IQ, Cu_UHQ, deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(Cu_IQ, Cu_UHQ,'kx', label='Messwerte')
x_plot = np.linspace(0, 10, 1000)
plt.plot(x_plot, x_plot*par[0]+par[1], 'r-', label="Lineare Regression")
plt.yticks([0, 5*1e-06, 10*1e-06, 15*1e-06, 20*1e-06],
           [0, 5, 10, 15, 20])
plt.legend(loc="best")
plt.xlabel(r'Quellstrom $I_Q\:/\:A$')
plt.ylabel(r'Hallspannung $U_H\:/\:\mu V$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotCu_IQ.pdf')
plt.close()

parCu_IQ=unp.uarray(par, err)

#Plot von Silber mit konstantem Querstrom %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B= (parhyauf[0].n + parhyab[0].n)/2*Ag_IB

par, covm = np.polyfit(B, Ag_UHB, deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(B, Ag_UHB,'kx', label='Messwerte')
x_plot = np.linspace(0, 1.3, 1000)
plt.plot(x_plot, x_plot*par[0]+par[1], 'r-', label="Lineare Regression")
plt.yticks([-170*1e-06, -175*1e-06, -180*1e-06, -185*1e-06, -190*1e-06, -195*1e-06],
           [-170, -175, -180, -185, -190, -195])
plt.legend(loc="best")
plt.xlabel(r'magnetische Feldstärke $B\:/\:T$')
plt.ylabel(r'Hallspannung $U_H\:/\:\mu V$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotAg_IB.pdf')
plt.close()

parAg_IB=unp.uarray(par, err)

#Plot von Silber mit konstantem B-Feld   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par, covm = np.polyfit(Ag_IQ, Ag_UHQ, deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(Ag_IQ, Ag_UHQ,'kx', label='Messwerte')
x_plot = np.linspace(0, 10, 1000)
plt.plot(x_plot, x_plot*par[0]+par[1], 'r-', label="Lineare Regression")
plt.yticks([0, -25*1e-06, -50*1e-06, -75*1e-06, -100*1e-06, -125*1e-06, -150*1e-06, -175*1e-06, -200*1e-06],
           [0, -25, -50, -75, -100, -125, -150, -175, -200])
plt.legend(loc="best")
plt.xlabel(r'Quellstrom $I_Q\:/\:A$')
plt.ylabel(r'Hallspannung $U_H\:/\:\mu V$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotAg_IQ.pdf')
plt.close()

parAg_IQ=unp.uarray(par, err)

#Plot von Zink mit konstantem Querstrom %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B= (parhyauf[0].n + parhyab[0].n)/2*Zn_IB

par, covm = np.polyfit(B, Zn_UHB, deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(B, Zn_UHB,'kx', label='Messwerte')
x_plot = np.linspace(0, 1.3, 1000)
plt.plot(x_plot, x_plot*par[0]+par[1], 'r-', label="Lineare Regression")
plt.yticks([ -313*1e-06, -315*1e-06, -317*1e-06, -319*1e-06],
           [ -313, -315, -317, -319])
plt.legend(loc="best")
plt.xlabel(r'magnetische Feldstärke $B\:/\:T$')
plt.ylabel(r'Hallspannung $U_H\:/\:\mu V$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotZn_IB.pdf')
plt.close()

parZn_IB=unp.uarray(par, err)

#Plot von Zink mit konstantem B-Feld %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

par, covm = np.polyfit(Zn_IQ, Zn_UHQ, deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(Zn_IQ, Zn_UHQ,'kx', label='Messwerte')
x_plot = np.linspace(0, 9, 1000)
plt.plot(x_plot, x_plot*par[0]+par[1], 'r-', label="Lineare Regression")
plt.yticks([0, -50*1e-06, -100*1e-06, -150*1e-06, -200*1e-06, -250*1e-06, -300*1e-06, -350*1e-06],
           [0, -50, -100, -150, -200, -250, -300, -350])
plt.legend(loc="best")
plt.xlabel(r'Quellstrom $I_Q\:/\:A$')
plt.ylabel(r'Hallspannung $U_H\:/\:\mu V$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotZn_IQ.pdf')
plt.close()

parZn_IQ=unp.uarray(par, err)

# Berechnungen relevanter Größen ################################################################################################

# Für Kupfer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# spezifischer Widerstand 

rho_Cu = rho(R_Cu, l_Cu, np.pi*(d_Cu_Draht/2)**2)

# Methode bei konstanten Querstrom /////////////////////////////////////////////////////////////////////////////

# Ladungsträger pro Volumen

n_Cu_IB = n(parCu_IB[0], d_Cu_Folie, 10)

# Hallkonstante

AH_Cu_IB = AH(n_Cu_IB)
print(AH_Cu_IB)
# Zahl der Ladungsträger pro Atom 

z_Cu_IB = z(n_Cu_IB, varrho_Cu, atom_mass_Cu)
print(z_Cu_IB)
# mittlere Flugzeit

tau_Cu_IB = tau(n_Cu_IB, varrho_Cu)

# mittlere Driftgeschwindigkeit

v_d_Cu_IB = v_d(n_Cu_IB)

# Beweglichkeit

mu_Cu_IB = mu(v_d_Cu_IB, n_Cu_IB, tau_Cu_IB)

# Totalgeschwindigkeit

v_Cu_IB = v(n_Cu_IB)

# mittlere freie Weglänge

l_Cu_IB = l(tau_Cu_IB, v_Cu_IB)

# Methode bei konstantem B-Feld //////////////////////////////////////////////////////////////////////////////////

# Ladungsträger pro Volumen

n_Cu_IQ = n(parCu_IQ[0], d_Cu_Folie, B1(parhyauf[0], parhyab[0],  5))

# Hallkonstante

AH_Cu_IQ = AH(n_Cu_IQ)
print(AH_Cu_IQ)
# Zahl der Ladungsträger pro Atom 

z_Cu_IQ = z(n_Cu_IQ, varrho_Cu, atom_mass_Cu)
print(z_Cu_IQ)
# mittlere Flugzeit

tau_Cu_IQ = tau(n_Cu_IQ, varrho_Cu)

# mittlere Driftgeschwindigkeit

v_d_Cu_IQ = v_d(n_Cu_IQ)

# Beweglichkeit

mu_Cu_IQ = mu(v_d_Cu_IQ, n_Cu_IQ, tau_Cu_IQ)

# Totalgeschwindigkeit

v_Cu_IQ = v(n_Cu_IQ)

# mittlere freie Weglänge

l_Cu_IQ = l(tau_Cu_IQ, v_Cu_IQ)

# Für Silber %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# spezifischer Widerstand 

rho_Ag = rho(R_Ag, l_Ag, np.pi*(d_Ag_Draht/2)**2)

# Methode bei konstanten Querstrom ///////////////////////////////////////////////////////////////////////////////////////

# Ladungsträger pro Volumen

n_Ag_IB = n(parAg_IB[0], d_Ag_Folie, 10)

# Hallkonstante

AH_Ag_IB = AH(n_Ag_IB)
print(AH_Ag_IB)
# Zahl der Ladungsträger pro Atom 

z_Ag_IB = z(n_Ag_IB, varrho_Ag, atom_mass_Ag)
print(z_Ag_IB)
# mittlere Flugzeit

tau_Ag_IB = tau(n_Ag_IB, varrho_Ag)

# mittlere Driftgeschwindigkeit

v_d_Ag_IB = v_d(n_Ag_IB)

# Beweglichkeit

mu_Ag_IB = mu(v_d_Ag_IB, n_Ag_IB, tau_Ag_IB)

# Totalgeschwindigkeit

v_Ag_IB = v(n_Ag_IB)

# mittlere freie Weglänge

l_Ag_IB = l(tau_Ag_IB, v_Ag_IB)

# Methode bei konstantem B-Feld //////////////////////////////////////////////////////////////////////////////////////

# Ladungsträger pro Volumen

n_Ag_IQ = n(parAg_IQ[0], d_Ag_Folie, B1(parhyauf[0], parhyab[0],  5))

# Hallkonstante

AH_Ag_IQ = AH(n_Ag_IQ)
print(AH_Ag_IQ)
# Zahl der Ladungsträger pro Atom

z_Ag_IQ = z(n_Ag_IQ, varrho_Ag, atom_mass_Ag)
print(z_Ag_IQ)
# mittlere Flugzeit

tau_Ag_IQ = tau(n_Ag_IQ, varrho_Ag)

# mittlere Driftgeschwindigkeit

v_d_Ag_IQ = v_d(n_Ag_IQ)

# Beweglichkeit

mu_Ag_IQ = mu(v_d_Ag_IQ, n_Ag_IQ, tau_Ag_IQ)

# Totalgeschwindigkeit

v_Ag_IQ = v(n_Ag_IQ)

# mittlere freie Weglänge

l_Ag_IQ = l(tau_Ag_IQ, v_Ag_IQ)

# Für Zink %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Methode bei konstanten Querstrom ///////////////////////////////////////////////////////////////////////

# Ladungsträger pro Volumen

n_Zn_IB = n(parZn_IB[0], d_Zn_Folie, 8)

# Hallkonstante

AH_Zn_IB = AH(n_Zn_IB)
print(AH_Zn_IB)
# Zahl der Ladungsträger pro Atom 

z_Zn_IB = z(n_Zn_IB, varrho_Zn, atom_mass_Zn)
print(z_Zn_IB)
# mittlere Flugzeit

tau_Zn_IB = tau(n_Zn_IB, varrho_Zn)

# mittlere Driftgeschwindigkeit

v_d_Zn_IB = v_d(n_Zn_IB)

# Beweglichkeit

mu_Zn_IB = mu(v_d_Zn_IB, n_Zn_IB, tau_Zn_IB)

# Totalgeschwindigkeit

v_Zn_IB = v(n_Zn_IB)

# mittlere freie Weglänge

l_Zn_IB = l(tau_Zn_IB, v_Zn_IB)

# Methode bei konstantem B-Feld //////////////////////////////////////////////////////////////////////////////////////

# Ladungsträger pro Volumen

n_Zn_IQ = n(parZn_IQ[0], d_Zn_Folie, B1(parhyauf[0], parhyab[0],  5)) 

# Hallkonstante

AH_Zn_IQ = AH(n_Zn_IQ)
print(AH_Zn_IQ)
# Zahl der Ladungsträger pro Atom 

z_Zn_IQ = z(n_Zn_IQ, varrho_Ag, atom_mass_Ag)
print(z_Zn_IQ)
# mittlere Flugzeit

tau_Zn_IQ = tau(n_Zn_IQ, varrho_Zn)

# mittlere Driftgeschwindigkeit

v_d_Zn_IQ = v_d(n_Zn_IB)

# Beweglichkeit

mu_Zn_IQ = mu(v_d_Zn_IQ, n_Zn_IQ, tau_Zn_IQ)

# Totalgeschwindigkeit

v_Zn_IQ = v(n_Zn_IQ)

# mittlere freie Weglänge

l_Zn_IQ = l(tau_Zn_IQ, v_Zn_IQ)

# Tex ########################################################################################################################################################



#
#
# Es fehlen noch die richtigen Größeneinheiten/Exponenten !!!
#
#



# Kupfer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# tex file for rho_Cu.tex

with open('build/rho_Cu.tex', 'w') as f:
  f.write(make_SI(rho_Cu,r'', figures=2))

# Konstanter Querstrom ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# tex file for n_Cu_IB.tex

with open('build/n_Cu_IB.tex', 'w') as f:
  f.write(make_SI(n_Cu_IB,r'', figures=2))
  
# tex file for AH_Cu_IB.tex

with open('build/AH_Cu_IB.tex', 'w') as f:
  f.write(make_SI(AH_Cu_IB,r'', figures=2))

# tex file for z_Cu_IB.tex

with open('build/z_Cu_IB.tex', 'w') as f:
  f.write(make_SI(z_Cu_IB,r'', figures=2))

# tex file for tau_Cu_IB.tex

with open('build/tau_Cu_IB.tex', 'w') as f:
  f.write(make_SI(tau_Cu_IB,r'', figures=2))

# tex file for v_d_Cu_IB.tex

with open('build/v_d_Cu_IB.tex', 'w') as f:
  f.write(make_SI(v_d_Cu_IB,r'', figures=2))

# tex file for mu_Cu_IB.tex

with open('build/mu_Cu_IB.tex', 'w') as f:
  f.write(make_SI(mu_Cu_IB,r'', figures=2))

# tex file for v_Cu_IB.tex

with open('build/v_Cu_IB.tex', 'w') as f:
  f.write(make_SI(v_Cu_IB,r'', figures=2))

# tex file for l_Cu_IB.tex

with open('build/l_Cu_IB.tex', 'w') as f:
  f.write(make_SI(l_Cu_IB,r'', figures=2))

# Konstantes B-Feld ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# tex file for n_Cu_IQ.tex

with open('build/n_Cu_IQ.tex', 'w') as f:
  f.write(make_SI(n_Cu_IQ,r'', figures=2))
  
# tex file for AH_Cu_IQ.tex

with open('build/AH_Cu_IQ.tex', 'w') as f:
  f.write(make_SI(AH_Cu_IQ,r'', figures=2))

# tex file for z_Cu_IQ.tex

with open('build/z_Cu_IQ.tex', 'w') as f:
  f.write(make_SI(z_Cu_IQ,r'', figures=2))

# tex file for tau_Cu_IQ.tex

with open('build/tau_Cu_IQ.tex', 'w') as f:
  f.write(make_SI(tau_Cu_IQ,r'', figures=2))

# tex file for v_d_Cu_IQ.tex

with open('build/v_d_Cu_IQ.tex', 'w') as f:
  f.write(make_SI(v_d_Cu_IQ,r'', figures=2))

# tex file for mu_Cu_IQ.tex

with open('build/mu_Cu_IQ.tex', 'w') as f:
  f.write(make_SI(mu_Cu_IQ,r'', figures=2))

# tex file for v_Cu_IQ.tex

with open('build/v_Cu_IQ.tex', 'w') as f:
  f.write(make_SI(v_Cu_IQ,r'', figures=2))

# tex file for l_Cu_IQ.tex

with open('build/l_Cu_IQ.tex', 'w') as f:
  f.write(make_SI(l_Cu_IQ,r'', figures=2))

# Silber %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# tex file for rho_Ag.tex

with open('build/rho_Ag.tex', 'w') as f:
  f.write(make_SI(rho_Ag,r'', figures=2))

# Konstanter Querstrom ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# tex file for n_Ag_IB.tex

with open('build/n_Ag_IB.tex', 'w') as f:
  f.write(make_SI(n_Ag_IB,r'', figures=2))
  
# tex file for AH_Ag_IB.tex

with open('build/AH_Ag_IB.tex', 'w') as f:
  f.write(make_SI(AH_Ag_IB,r'', figures=2))

# tex file for z_Ag_IB.tex

with open('build/z_Ag_IB.tex', 'w') as f:
  f.write(make_SI(z_Ag_IB,r'', figures=2))

# tex file for tau_Ag_IB.tex

with open('build/tau_Ag_IB.tex', 'w') as f:
  f.write(make_SI(tau_Ag_IB,r'', figures=2))

# tex file for v_d_Ag_IB.tex

with open('build/v_d_Ag_IB.tex', 'w') as f:
  f.write(make_SI(v_d_Ag_IB,r'', figures=2))

# tex file for mu_Ag_IB.tex

with open('build/mu_Ag_IB.tex', 'w') as f:
  f.write(make_SI(mu_Ag_IB,r'', figures=2))

# tex file for v_Ag_IB.tex

with open('build/v_Ag_IB.tex', 'w') as f:
  f.write(make_SI(v_Ag_IB,r'', figures=2))

# tex file for l_Ag_IB.tex

with open('build/l_Ag_IB.tex', 'w') as f:
  f.write(make_SI(l_Ag_IB,r'', figures=2))

# Konstantes B-Feld ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# tex file for n_Ag_IQ.tex

with open('build/n_Ag_IQ.tex', 'w') as f:
  f.write(make_SI(n_Ag_IQ,r'', figures=2))
  
# tex file for AH_Ag_IQ.tex

with open('build/AH_Ag_IQ.tex', 'w') as f:
  f.write(make_SI(AH_Ag_IQ,r'', figures=2))

# tex file for z_Ag_IQ.tex

with open('build/z_Ag_IQ.tex', 'w') as f:
  f.write(make_SI(z_Ag_IQ,r'', figures=2))

# tex file for tau_Ag_IQ.tex

with open('build/tau_Ag_IQ.tex', 'w') as f:
  f.write(make_SI(tau_Ag_IQ,r'', figures=2))

# tex file for v_d_Ag_IQ.tex

with open('build/v_d_Ag_IQ.tex', 'w') as f:
  f.write(make_SI(v_d_Ag_IQ,r'', figures=2))

# tex file for mu_Ag_IQ.tex

with open('build/mu_Ag_IQ.tex', 'w') as f:
  f.write(make_SI(mu_Ag_IQ,r'', figures=2))

# tex file for v_Ag_IQ.tex

with open('build/v_Ag_IQ.tex', 'w') as f:
  f.write(make_SI(v_Ag_IQ,r'', figures=2))

# tex file for l_Ag_IQ.tex

with open('build/l_Ag_IQ.tex', 'w') as f:
  f.write(make_SI(l_Ag_IQ,r'', figures=2))

# Zink %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Konstanter Querstrom ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# tex file for n_Zn_IB.tex

with open('build/n_Zn_IB.tex', 'w') as f:
  f.write(make_SI(n_Zn_IB,r'', figures=2))
  
# tex file for AH_Zn_IB.tex

with open('build/AH_Zn_IB.tex', 'w') as f:
  f.write(make_SI(AH_Zn_IB,r'', figures=2))

# tex file for z_Zn_IB.tex

with open('build/z_Zn_IB.tex', 'w') as f:
  f.write(make_SI(z_Zn_IB,r'', figures=2))

# tex file for tau_Zn_IB.tex

with open('build/tau_Zn_IB.tex', 'w') as f:
  f.write(make_SI(tau_Zn_IB,r'', figures=2))

# tex file for v_d_Zn_IB.tex

with open('build/v_d_Zn_IB.tex', 'w') as f:
  f.write(make_SI(v_d_Zn_IB,r'', figures=2))

# tex file for mu_Zn_IB.tex

with open('build/mu_Zn_IB.tex', 'w') as f:
  f.write(make_SI(mu_Zn_IB,r'', figures=2))

# tex file for v_Ag_IB.tex

with open('build/v_Ag_IB.tex', 'w') as f:
  f.write(make_SI(v_Ag_IB,r'', figures=2))

# tex file for l_Zn_IB.tex

with open('build/l_Zn_IB.tex', 'w') as f:
  f.write(make_SI(l_Zn_IB,r'', figures=2))

# Konstantes B-Feld ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# tex file for n_Zn_IQ.tex

with open('build/n_Zn_IQ.tex', 'w') as f:
  f.write(make_SI(n_Zn_IQ,r'', figures=2))
  
# tex file for AH_Zn_IQ.tex

with open('build/AH_Zn_IQ.tex', 'w') as f:
  f.write(make_SI(AH_Zn_IQ,r'', figures=2))

# tex file for z_Zn_IQ.tex

with open('build/z_Zn_IQ.tex', 'w') as f:
  f.write(make_SI(z_Zn_IQ,r'', figures=2))

# tex file for tau_Zn_IQ.tex

with open('build/tau_Zn_IQ.tex', 'w') as f:
  f.write(make_SI(tau_Zn_IQ,r'', figures=2))

# tex file for v_d_Zn_IQ.tex

with open('build/v_d_Zn_IQ.tex', 'w') as f:
  f.write(make_SI(v_d_Zn_IQ,r'', figures=2))

# tex file for mu_Zn_IQ.tex

with open('build/mu_Zn_IQ.tex', 'w') as f:
  f.write(make_SI(mu_Zn_IQ,r'', figures=2))

# tex file for v_Zn_IQ.tex

with open('build/v_Zn_IQ.tex', 'w') as f:
  f.write(make_SI(v_Zn_IQ,r'', figures=2))

# tex file for l_Zn_IQ.tex

with open('build/l_Zn_IQ.tex', 'w') as f:
  f.write(make_SI(l_Zn_IQ,r'', figures=2))

# Tabellen ########################################################################################################################################################

# Kupfer ===============================================================================================================================================================

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    \multicolumn{2}{c}{$\text{Konstanter Querstrom =} \; \SI{10}{\ampere}$} & \multicolumn{2}{c}{$\text{Konstantes B-Feld =} \; \SI{5}{\ampere}$}\\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4}
    \multicolumn{1}{c}{$\text{Strom} \; I_B \:/\: \si{\ampere}$} & \multicolumn{1}{c}{$\text{Spannung} \; U_H \:/\: \si{\milli\volt} $} & 
    \multicolumn{1}{c}{$\text{Strom} \; I_Q \:/\: \si{\ampere}$} & \multicolumn{1}{c}{$\text{Spannung} \; U_H \:/\: \si{\milli\volt}$}\\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4} 

'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.2f} & {2:1.0f} & {3:1.2f}\\'

with open('build/Cu_table.tex', 'w') as g:
    g.write(table_header)
    for row in zip(Cu_IB, Cu_UHB*1e+04, Cu_IQ, Cu_UHQ*1e+04):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

# Zink ===============================================================================================================================================================

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    \multicolumn{2}{c}{$\text{Konstanter Querstrom =} \; \SI{8}{\ampere}$} & \multicolumn{2}{c}{$\text{Konstantes B-Feld =} \; \SI{5}{\ampere}$}\\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4}
    \multicolumn{1}{c}{$\text{Strom} \; I_B \:/\: \si{\ampere}$} & \multicolumn{1}{c}{$\text{Spannung} \; U_H \:/\: \si{\milli\volt} $} & 
    \multicolumn{1}{c}{$\text{Strom} \; I_Q \:/\: \si{\ampere}$} & \multicolumn{1}{c}{$\text{Spannung} \; U_H \:/\: \si{\milli\volt}$}\\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4} 

'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.2f} & {2:1.1f} & {3:1.2f}\\'

with open('build/Zn_table.tex', 'w') as g:
    g.write(table_header)
    for row in zip(Zn_IB, Zn_UHB*1e+04, Zn_IQ, Zn_UHQ*1e+04):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

# Silber ===============================================================================================================================================================

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    \multicolumn{2}{c}{$\text{Konstanter Querstrom =} \; \SI{10}{\ampere}$} & \multicolumn{2}{c}{$\text{Konstantes B-Feld =} \; \SI{5}{\ampere}$}\\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4}
    \multicolumn{1}{c}{$\text{Strom} \; I_B \:/\: \si{\ampere}$} & \multicolumn{1}{c}{$\text{Spannung} \; U_H \:/\: \si{\milli\volt} $} & 
    \multicolumn{1}{c}{$\text{Strom} \; I_Q \:/\: \si{\ampere}$} & \multicolumn{1}{c}{$\text{Spannung} \; U_H \:/\: \si{\milli\volt}$}\\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4} 

'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.2f} & {2:1.0f} & {3:1.2f}\\'

with open('build/Ag_table.tex', 'w') as g:
    g.write(table_header)
    for row in zip(Ag_IB, Ag_UHB*1e+04, Ag_IQ, Ag_UHQ*1e+04):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

# Hysterese ===============================================================================================================================================================

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    \multicolumn{2}{c}{$\text{Aufsteigend von} \; \SI{0}{\ampere}$} & \multicolumn{2}{c}{$\text{Absteigend von} \; \SI{4.5}{\ampere}$}\\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4}
    \multicolumn{1}{c}{$\text{Strom} \; I_{auf} \:/\: \si{\ampere}$} & \multicolumn{1}{c}{$\text{Magnetfeld} \; B_{auf} \:/\: \si{\milli\tesla} $} & 
    \multicolumn{1}{c}{$\text{Strom} \; I_{ab} \:/\: \si{\ampere}$} & \multicolumn{1}{c}{$\text{Magnetfeld} \; B_{ab} \:/\: \si{\milli\tesla}$}\\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4} 

'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.1f} & {2:1.1f} & {3:1.1f}\\'

with open('build/Hysterese_table.tex', 'w') as g:
    g.write(table_header)
    for row in zip(hy_Iauf, hy_Bauf*1e+03, hy_Iab, hy_Bab*1e+03):
        g.write(row_template.format(*row).replace('nan',' ') )
        g.write('\n')
    g.write(table_footer)