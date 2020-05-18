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


def rel_err(mess, lit):
    return abs((unp.nominal_values(mess)-unp.nominal_values(lit))/unp.nominal_values(lit))*100

def gauß(grad, A, B, C, D):
    return A* np.exp(-((grad-B)/C)**2)+D

def E(grad):
  return const.h*const.c/(2*201.4e-12*unp.sin(grad)*const.e)

def I(max, min):
  return min + (max-min)/2

def sigma_1(E_abs, z):
  return  z-unp.sqrt(E_abs/(13.6))

def sigma_2(E ,E_abs ,z ,ml):
  return z-unp.sqrt(E_abs/(13.6)-E/(13.6))*ml

alpha = 1/(137.302)

def sigma(E, z):
  return z- unp.sqrt(E/13.6-alpha**2*z**4/4)





# Messwerte ##########################################################################################################################


# Messwerte für die Überprüfung der Bragg-Bedingung

theta_deg_Bragg, N_Bragg = np.genfromtxt('Bragg.dat', unpack=True)


# Messwerte für das Absorptionsspektrum von Brom

theta_deg_Br, N_Br = np.genfromtxt('Brom.dat', unpack=True)


# Messwerte für das Emmissionsspektrum von Kupfer

theta_deg_Cu, N_Cu = np.genfromtxt('Emissionsspektrum.dat', unpack=True)


# Messwerte für das Absorptionsspektrum von Gallium

theta_deg_Ga, N_Ga = np.genfromtxt('Gallium.dat', unpack=True)


# Messwerte für das Absorptionsspektrum von Rubidium

theta_deg_Rb, N_Rb = np.genfromtxt('Rubidium.dat', unpack=True)


# Messwerte für das Absorptionsspektrum von Strontium

theta_deg_Sr, N_Sr = np.genfromtxt('Strontium.dat', unpack=True)


# Messwerte für das Absorptionsspektrum von Zink

theta_deg_Zn, N_Zn = np.genfromtxt('Zink.dat', unpack=True)


# Messwerte für das Absorptionsspektrum von Zirkonium

theta_deg_Zr, N_Zr = np.genfromtxt('Zirkonium.dat', unpack=True)



# Plots ################################################################################################################


# Plot zu Überprüfung der Bragg-Bedingung

popBragg, pcovBragg = curve_fit(
    gauß,
    theta_deg_Bragg[:41],
    N_Bragg[:41],
    sigma=None,
    absolute_sigma=True,
    p0=[150,28,-30, 20]
    )

err = np.sqrt(np.diag(pcovBragg))

plt.plot(theta_deg_Bragg[:41], N_Bragg[:41],'kx', label='Messwerte')
x_plot = np.linspace(26, 30, 1000)
plt.plot(x_plot, gauß(x_plot, *popBragg), 'r-', label='Fitkurve')
plt.axvline(x=28.15, color='r', linestyle=':', label='$Maximum$')
plt.legend(loc="best")
plt.xlabel(r'Winkel $\theta \:/\:°$')
plt.ylabel(r'Intensität $I\:/\:Imp/s$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_Bragg.pdf')
plt.close()

par = unp.uarray(popBragg, err)

# Plot des Kupfer Emissionsspektrums

plt.plot(theta_deg_Cu, N_Cu,'kx', label='Messwerte')
plt.axvline(x=20.2, color='r', linestyle=':', label='$K_{\\beta}- Linie$')
plt.axvline(x=22.5, color='m', linestyle=':', label='$K_{\\alpha}-Linie$')
plt.axhline(y=1599/2 ,color='g', linestyle='-', label='$\Delta E_{FWHM}\; K_{\\beta} $' )
plt.axhline(y=2525 ,color='b', linestyle='-', label='$\Delta E_{FWHM}\; K_{\\alpha} $' )
plt.legend(loc="best")
plt.annotate('$K_{\\beta}$', xy=(20.2, 1700), size = 15

)
plt.annotate('$K_{\\alpha} $', xy=(23, 5050), size = 15

)
plt.annotate('Bremsspektrum', xy=(10, 800), size = 20

)
plt.xlabel(r'Winkel $\theta \:/\:°$')
plt.ylabel(r'Intensität $I\:/\:Imp/s$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_Cu.pdf')
plt.close()


# Plot zu Brom-Absorptionsspektrum

plt.plot(theta_deg_Br, N_Br,'kx', label='Messwerte')
plt.axhline(y=I(27, 9), color='b', linestyle='-', label='Mitte der K-Kante' )
plt.axvline(x=13.2, color='r', linestyle='-', label='Winkel der K-Kante')
plt.legend(loc="best")
plt.xlabel(r'Winkel $\theta \:/\:°$')
plt.ylabel(r'Intensität $I\:/\:Imp/s$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_Br.pdf')
plt.close()


# Plot zu Gallium-Absorptionsspektrum

plt.plot(theta_deg_Ga, N_Ga,'kx', label='Messwerte')
plt.axhline(y=I(122, 66), color='b', linestyle='-', label='Mitte der K-Kante' )
plt.axvline(x=17.35, color='r', linestyle='-', label='Winkel der K-Kante')
plt.legend(loc="best")
plt.xlabel(r'Winkel $\theta \:/\:°$')
plt.ylabel(r'Intensität $I\:/\:Imp/s$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_Ga.pdf')
plt.close()


# Plot zu Rubidium-Absorptionsspektrum

plt.plot(theta_deg_Rb, N_Rb,'kx', label='Messwerte')
plt.axhline(y=I(64, 10), color='b', linestyle='-', label='Mitte der K-Kante' )
plt.axvline(x=11.77, color='r', linestyle='-', label='Winkel der K-Kante')
plt.legend(loc="best")
plt.xlabel(r'Winkel $\theta \:/\:°$')
plt.ylabel(r'Intensität $I\:/\:Imp/s$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_Rb.pdf')
plt.close()


# Plot zu Strontium-Absorptionsspektrum

plt.plot(theta_deg_Sr, N_Sr,'kx', label='Messwerte')
plt.axhline(y=I(193, 40), color='b', linestyle='-', label='Mitte der K-Kante' )
plt.axvline(x=11.08, color='r', linestyle='-', label='Winkel der K-Kante')
plt.legend(loc="best")
plt.xlabel(r'Winkel $\theta \:/\:°$')
plt.ylabel(r'Intensität $I\:/\:Imp/s$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_Sr.pdf')
plt.close()


# Plot zu Zink-Absorptionsspektrum

plt.plot(theta_deg_Zn, N_Zn,'kx', label='Messwerte')
plt.axhline(y=I(102, 54), color='b', linestyle='-', label='Mitte der K-Kante' )
plt.axvline(x=18.65, color='r', linestyle='-', label='Winkel der K-Kante')
plt.legend(loc="best")
plt.xlabel(r'Winkel $\theta \:/\:°$')
plt.ylabel(r'Intensität $I\:/\:Imp/s$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_Zn.pdf')
plt.close()


# Plot zu Zirkonium-Absorptionsspektrum

plt.plot(theta_deg_Zr, N_Zr,'kx', label='Messwerte')
plt.axhline(y=I(301, 112), color='b', linestyle='-', label='Mitte der K-Kante' )
plt.axvline(x=9.95, color='r', linestyle='-', label='Winkel der K-Kante')
plt.legend(loc="best")
plt.xlabel(r'Winkel $\theta \:/\:°$')
plt.ylabel(r'Intensität $I\:/\:Imp/s$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_Zr.pdf')
plt.close()

# Tex-Files ########################################################################################################################

#tex file of A_Bragg  ###########################################################################

with open('build/A_Bragg.tex', 'w') as f: 
  f.write(make_SI(par[0] ,r'\hertz' ,figures=1))

with open('build/B_Bragg.tex', 'w') as f: 
  f.write(make_SI(par[1] ,r'' ,figures=1))

with open('build/C_Bragg.tex', 'w') as f: 
  f.write(make_SI(par[2] ,r'' ,figures=1))

with open('build/D_Bragg.tex', 'w') as f: 
  f.write(make_SI(par[3] ,r'\hertz' ,figures=1))

# tex files for Kupfer-Emissionsspektrum #########################################
theta_alpha =ufloat(np.deg2rad(22.5), np.deg2rad(0.1))
theta_beta =ufloat(np.deg2rad(20.2), np.deg2rad(0.1))

with open('build/E_alpha.tex', 'w') as f: 
  f.write(make_SI(E(theta_alpha) ,r'\electronvolt' ,figures=1))

with open('build/E_beta.tex', 'w') as f: 
  f.write(make_SI(E(theta_beta) ,r'\electronvolt' ,figures=1))

# tex file for Halbwertsbreite

with open('build/Delta_E_alpha.tex', 'w') as f: 
  f.write(make_SI(abs(E(ufloat(np.deg2rad(22.85),np.deg2rad(0.1)))-E(ufloat(np.deg2rad(22.35),np.deg2rad(0.1)))) ,r'\electronvolt' ,figures=1))

with open('build/Delta_E_beta.tex', 'w') as f: 
  f.write(make_SI(abs(E(ufloat(np.deg2rad(20.55),np.deg2rad(0.1)))-E(ufloat(np.deg2rad(20.05),np.deg2rad(0.1)))) ,r'\electronvolt' ,figures=1))


# tex file for Auflösung

with open('build/A_alpha.tex', 'w') as f: 
  f.write(make_SI(E(theta_alpha)/abs(E(ufloat(np.deg2rad(22.85),np.deg2rad(0.1)))-E(ufloat(np.deg2rad(22.35),np.deg2rad(0.1)))) ,r'' ,figures=1))

with open('build/A_beta.tex', 'w') as f: 
  f.write(make_SI(E(theta_beta)/abs(E(ufloat(np.deg2rad(20.55),np.deg2rad(0.1)))-E(ufloat(np.deg2rad(20.05),np.deg2rad(0.1)))) ,r'' ,figures=1))


# tex file for sigma _1

E_abs_Cu = 8987.96
z_Cu = 29


with open('build/sigma_1_Cu.tex', 'w') as f: 
  f.write(make_SI(sigma_1(E_abs_Cu ,z_Cu) ,r'' ,figures=2))


# tex file for sigma_2,3
E_alpha=E(theta_alpha)
E_beta=E(theta_beta)

with open('build/sigma_2.tex', 'w') as f: 
  f.write(make_SI(sigma_2(E_alpha ,E_abs_Cu ,z_Cu ,2) ,r'' ,figures=2))

sigma_alpha=sigma_2(E_alpha ,E_abs_Cu ,z_Cu ,2)

with open('build/sigma_3.tex', 'w') as f: 
  f.write(make_SI(sigma_2(E_beta ,E_abs_Cu ,z_Cu ,3) ,r'' ,figures=2))

sigma_beta=sigma_2(E_beta ,E_abs_Cu ,z_Cu ,3)

# tex file for sigma_2,3 theo

with open('build/sigma_2_theo.tex', 'w') as f: 
  f.write(make_SI(sigma_2(8048.11 ,E_abs_Cu ,z_Cu ,2) ,r'' ,figures=2))


with open('build/sigma_3_theo.tex', 'w') as f: 
  f.write(make_SI(sigma_2(8906.9 ,E_abs_Cu ,z_Cu ,3) ,r'' ,figures=2))


# tex file for Brom

I_Brom=I(27, 9)


with open('build/I_Brom.tex', 'w') as f: 
  f.write(make_SI(I_Brom ,r'\hertz' ,figures=2))


Brom_rad = ufloat(np.deg2rad(13.2),np.deg2rad(0.1))

E_Brom=E(Brom_rad)

with open('build/E_Brom.tex', 'w') as f: 
  f.write(make_SI(E_Brom ,r'\electronvolt' ,figures=2))

sigma_Brom=sigma(E_Brom, 35)

with open('build/sigma_Brom.tex', 'w') as f: 
  f.write(make_SI(sigma(E_Brom, 35) ,r'' ,figures=2))

E_Brom_theo =13483.86

sigma_Brom_theo=sigma(E_Brom_theo, 35)
 
with open('build/sigma_Brom_theo.tex', 'w') as f: 
  f.write(make_SI(sigma(E_Brom_theo, 35) ,r'' ,figures=2))


# tex file for Gallium

I_Gallium=I(122, 66)


with open('build/I_Gallium.tex', 'w') as f: 
  f.write(make_SI(I_Gallium ,r'\hertz' ,figures=2))


Gallium_rad = ufloat(np.deg2rad(17.35),np.deg2rad(0.1))

E_Gallium=E(Gallium_rad)

with open('build/E_Gallium.tex', 'w') as f: 
  f.write(make_SI(E_Gallium ,r'\electronvolt' ,figures=2))

sigma_Gallium=sigma(E_Gallium, 31)
with open('build/sigma_Gallium.tex', 'w') as f: 
  f.write(make_SI(sigma(E_Gallium, 31) ,r'' ,figures=2))

E_Gallium_theo =10377.76

sigma_Gallium_theo=sigma(E_Gallium_theo, 31)

with open('build/sigma_Gallium_theo.tex', 'w') as f: 
  f.write(make_SI(sigma(E_Gallium_theo, 31) ,r'' ,figures=2))


# tex file for Rubidium

I_Rb=I(64, 10)


with open('build/I_Rb.tex', 'w') as f: 
  f.write(make_SI(I_Rb ,r'\hertz' ,figures=2))


Rb_rad = ufloat(np.deg2rad(11.77),np.deg2rad(0.1))

E_Rb=E(Rb_rad)

with open('build/E_Rb.tex', 'w') as f: 
  f.write(make_SI(E_Rb ,r'\electronvolt' ,figures=2))

sigma_Rb=sigma(E_Rb, 37)
with open('build/sigma_Rb.tex', 'w') as f: 
  f.write(make_SI(sigma(E_Rb, 37) ,r'' ,figures=2))

E_Rb_theo =15207.74
sigma_Rb_theo=sigma(E_Rb_theo, 37)

with open('build/sigma_Rb_theo.tex', 'w') as f: 
  f.write(make_SI(sigma(E_Rb_theo, 37) ,r'' ,figures=2))


# tex file for Strontium

I_Sr=I(193, 40)


with open('build/I_Sr.tex', 'w') as f: 
  f.write(make_SI(I_Sr ,r'\hertz' ,figures=2))


Sr_rad = ufloat(np.deg2rad(11.08),np.deg2rad(0.1))

E_Sr=E(Sr_rad)

with open('build/E_Sr.tex', 'w') as f: 
  f.write(make_SI(E_Sr ,r'\electronvolt' ,figures=2))

sigma_Sr=sigma(E_Sr, 38)
with open('build/sigma_Sr.tex', 'w') as f: 
  f.write(make_SI(sigma(E_Sr, 38) ,r'' ,figures=2))

E_Sr_theo =16115.26
sigma_Sr_theo=sigma(E_Sr_theo, 38)

with open('build/sigma_Sr_theo.tex', 'w') as f: 
  f.write(make_SI(sigma(E_Sr_theo, 38) ,r'' ,figures=2))


# tex file for Zink

I_Zn=I(102, 54)


with open('build/I_Zn.tex', 'w') as f: 
  f.write(make_SI(I_Zn ,r'\hertz' ,figures=2))


Zn_rad = ufloat(np.deg2rad(18.65),np.deg2rad(0.1))

E_Zn=E(Zn_rad)

with open('build/E_Zn.tex', 'w') as f: 
  f.write(make_SI(E_Zn ,r'\electronvolt' ,figures=2))

sigma_Zn=sigma(E_Zn, 30)

with open('build/sigma_Zn.tex', 'w') as f: 
  f.write(make_SI(sigma(E_Zn, 30) ,r'' ,figures=2))

E_Zn_theo =9668.55

sigma_Zn_theo=sigma(E_Zn_theo, 30) 

with open('build/sigma_Zn_theo.tex', 'w') as f: 
  f.write(make_SI(sigma(E_Zn_theo, 30) ,r'' ,figures=2))


# tex file for Zirkonium

I_Zr=I(301, 112)


with open('build/I_Zr.tex', 'w') as f: 
  f.write(make_SI(I_Zr ,r'\hertz' ,figures=2))


Zr_rad = ufloat(np.deg2rad(9.95),np.deg2rad(0.1))

E_Zr=E(Zr_rad)

with open('build/E_Zr.tex', 'w') as f: 
  f.write(make_SI(E_Zr ,r'\electronvolt' ,figures=2))

sigma_Zr=sigma(E_Zr, 40)

with open('build/sigma_Zr.tex', 'w') as f: 
  f.write(make_SI(sigma(E_Zr, 40) ,r'' ,figures=2))

E_Zr_theo =18008.15

sigma_Zr_theo=sigma(E_Zr_theo, 40)


with open('build/sigma_Zr_theo.tex', 'w') as f: 
  f.write(make_SI(sigma(E_Zr_theo, 40) ,r'' ,figures=2))


z_eff=np.array([35-unp.nominal_values(sigma_Brom), 31-unp.nominal_values(sigma_Gallium), 37-unp.nominal_values(sigma_Rb), 38-unp.nominal_values(sigma_Sr), 30-unp.nominal_values(sigma_Zn), unp.nominal_values(40-sigma_Zr)])

z_eff_theo=np.array([35-unp.nominal_values(sigma_Brom_theo), 31-unp.nominal_values(sigma_Gallium_theo), 37-unp.nominal_values(sigma_Rb_theo), 38-unp.nominal_values(sigma_Sr_theo), 30-unp.nominal_values(sigma_Zn_theo), unp.nominal_values(40-sigma_Zr_theo)] )

E=np.array([unp.nominal_values(E_Brom), unp.nominal_values(E_Gallium), unp.nominal_values(E_Rb), unp.nominal_values(E_Sr), unp.nominal_values(E_Zn), unp.nominal_values(E_Zr)])

E_theo=np.array([unp.nominal_values(E_Brom_theo), unp.nominal_values(E_Gallium_theo), unp.nominal_values(E_Rb_theo), unp.nominal_values(E_Sr_theo), unp.nominal_values(E_Zn_theo), unp.nominal_values(E_Zr_theo)])


par, cov = np.polyfit(z_eff, np.sqrt(E), deg=1, cov =True)
err = np.sqrt(np.diag(cov))



plt.plot(z_eff,np.sqrt(E) ,'kx', label='Messwerte')
x_plot = np.linspace(26, 36, 10000)
plt.plot(x_plot, par[0]*x_plot+par[1], 'k-',label='Ausgleichsgerade' )
plt.plot(x_plot, np.sqrt(13.6)*x_plot, 'b-',label='Ausgleichsgerade Theorie' )
plt.legend(loc="best")
plt.xlabel(r'effektiven Kernladung $z_{eff}$ ')
plt.ylabel(r'Wurzel der Absorptionsenergie $\sqrt{E \:/\:eV}$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_Ez.pdf')
plt.close()

par = unp.uarray(par,err)


# tex file for parameters

with open('build/a.tex', 'w') as f: 
  f.write(make_SI(par[0] ,r'\electronvolt' ,figures=2))

with open('build/b.tex', 'w') as f: 
  f.write(make_SI(par[1] ,r'\electronvolt' ,figures=2))

with open('build/R.tex', 'w') as f: 
  f.write(make_SI(par[0]**2*const.e/const.h *1e-15,r'\peta\hertz' ,figures=2))

with open('build/R_theo.tex', 'w') as f: 
  f.write(make_SI(13.6*const.e/const.h *1e-15,r'\peta\hertz' ,figures=2))


#Relativen Fehler ##########################################################

# Maximum der Kurve

with open('build/relerr_max.tex', 'w') as f: 
  f.write(make_SI(rel_err(28.15,28),r'' ,figures=2))

# Emissionsspektrum Kupfer

with open('build/relerr_E_alpha.tex', 'w') as f: 
  f.write(make_SI(rel_err(E_alpha,8048.11),r'' ,figures=2))

with open('build/relerr_E_beta.tex', 'w') as f: 
  f.write(make_SI(rel_err(E_beta,8906.9),r'' ,figures=2))

with open('build/relerr_sigma_alpha.tex', 'w') as f: 
  f.write(make_SI(rel_err(sigma_alpha,12.37),r'' ,figures=2))


with open('build/relerr_sigma_beta.tex', 'w') as f: 
  f.write(make_SI(rel_err(sigma_beta,21.68),r'' ,figures=2))


# Absorptionsspektren

# Brom

with open('build/relerr_E_Brom.tex', 'w') as f: 
  f.write(make_SI(rel_err(E_Brom,E_Brom_theo),r'' ,figures=2))

with open('build/relerr_sigma_Brom.tex', 'w') as f: 
  f.write(make_SI(rel_err(sigma_Brom,sigma_Brom_theo),r'' ,figures=2))


# Gallium

with open('build/relerr_E_Gallium.tex', 'w') as f: 
  f.write(make_SI(rel_err(E_Gallium,E_Gallium_theo),r'' ,figures=2))

with open('build/relerr_sigma_Gallium.tex', 'w') as f: 
  f.write(make_SI(rel_err(sigma_Gallium,sigma_Gallium_theo),r'' ,figures=2))

# Rubidium

with open('build/relerr_E_Rb.tex', 'w') as f: 
  f.write(make_SI(rel_err(E_Rb,E_Rb_theo),r'' ,figures=2))

with open('build/relerr_sigma_Rb.tex', 'w') as f: 
  f.write(make_SI(rel_err(sigma_Rb,sigma_Rb_theo),r'' ,figures=2))
print(sigma_Rb)
print(sigma_Rb_theo)
# Strontium

with open('build/relerr_E_Sr.tex', 'w') as f: 
  f.write(make_SI(rel_err(E_Sr,E_Sr_theo),r'' ,figures=2))

with open('build/relerr_sigma_Sr.tex', 'w') as f: 
  f.write(make_SI(rel_err(sigma_Sr,sigma_Sr_theo),r'' ,figures=2))

# Zink

with open('build/relerr_E_Zn.tex', 'w') as f: 
  f.write(make_SI(rel_err(E_Zn,E_Zn_theo),r'' ,figures=2))

print(rel_err(E_Zn,E_Zn_theo))

with open('build/relerr_sigma_Zn.tex', 'w') as f: 
  f.write(make_SI(rel_err(sigma_Zn,sigma_Zn_theo),r'' ,figures=2))

# Zirkonium

with open('build/relerr_E_Zr.tex', 'w') as f: 
  f.write(make_SI(rel_err(E_Zr,E_Zr_theo),r'' ,figures=2))

with open('build/relerr_sigma_Zr.tex', 'w') as f: 
  f.write(make_SI(rel_err(sigma_Zr,sigma_Zr_theo),r'' ,figures=2))

# Rydberg-Frequenz


with open('build/relerr_R.tex', 'w') as f: 
  f.write(make_SI(rel_err(par[0]**2*const.e/const.h *1e-15, 13.6*const.e/const.h *1e-15),r'' ,figures=2))



# Tabellen #######################################################################################################################

# Tabelle für die Überprüfung der Bragg Bedingung

theta_deg_Bragg_1, theta_deg_Bragg_2 = np.array_split(theta_deg_Bragg, 2)
N_Bragg_1, N_Bragg_2 = np.array_split(N_Bragg, 2)

table_header = r'''
  \begin{tabular}{S[table-format=2.1] S[table-format=3.1] S[table-format=2.1] S[table-format=3.1]}
    \toprule
    \multicolumn{1}{c}{Winkel $\theta \:/\:° $ } & \multicolumn{1}{c}{Zählrate $N\:/\:Imp/s $ } & \multicolumn{1}{c}{Winkel $\theta \:/\:° $ } & \multicolumn{1}{c}{Zählrate $N\:/\:Imp/s $ } \\
    \cmidrule(lr){1-2} \cmidrule(lr){3-4}
'''
table_footer = r'''    \bottomrule

  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.1f} & {2:1.1f} & {3:1.1f}  \\'


with open('build/table_Bragg.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_deg_Bragg_1, N_Bragg_1, theta_deg_Bragg_2, N_Bragg_2):
        g.write(row_template.format(*row).replace('nan',' ').replace('.',','))
        g.write('\n')
    g.write(table_footer)



# Tabelle für das Emissionsspektrum von Kupfer

theta_deg_Cu_1, theta_deg_Cu_2, theta_deg_Cu_3 = np.array_split(theta_deg_Cu, 3)
N_Cu_1, N_Cu_2, N_Cu_3 = np.array_split(N_Cu, 3)


table_header = r'''
  \begin{longtable}{c c c c c c}
  \caption{Messwerte für das Emissionsspektrum von Kupfer}\\
    \toprule
    \multicolumn{1}{c}{Winkel $\theta \:/\:° $ } & \multicolumn{1}{c}{Zählrate $N\:/\:Imp/s $ } & \multicolumn{1}{c}{ $\theta \:/\:°$ } & \multicolumn{1}{c}{$N\:/\:Imp/s $ }& \multicolumn{1}{c}{ $\theta \:/\:° $}& \multicolumn{1}{c}{ $N\:/\:Imp/s $  }\\
    \cmidrule(lr){1-2} \cmidrule(lr{0,5em}){3-4} \cmidrule(lr{0,5em}){5-6}
'''
table_footer = r'''    \bottomrule
  \label{tab:2}
  \end{longtable}
'''
row_template = r'     {0:1.1f} & {1:1.0f} & {2:1.1f} & {3:1.0f} & {4:1.1f} & {5:1.0f} \\'


with open('build/table_Cu.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_deg_Cu_1, N_Cu_1, theta_deg_Cu_2, N_Cu_2, theta_deg_Cu_3, N_Cu_3):
        g.write(row_template.format(*row).replace('.',','))
        g.write('\n')
    g.write(table_footer)


# Tabelle für das Absorptionsspektrum für Brom

theta_deg_Br_1, theta_deg_Br_2 = np.array_split(theta_deg_Br, 2)
N_Br_1, N_Br_2 = np.array_split(N_Br, 2)


table_header = r'''
  \begin{longtable}{S[table-format=2.1] S[table-format=2.1] S[table-format=2.1] S[table-format=2.1]}
  \caption{Messwerte für das Absorptionsspektrum von Brom}\\
    \toprule
    \multicolumn{1}{c}{ $\theta \:/\:°$ } & \multicolumn{1}{c}{$N\:/\:Imp/s $ }& \multicolumn{1}{c}{ $\theta \:/\:° $}& \multicolumn{1}{c}{ $N\:/\:Imp/s $  }\\
    \cmidrule(lr){1-2} \cmidrule(lr{0,5em}){3-4}
'''
table_footer = r'''    \bottomrule
  \label{tab:3}
  \end{longtable}
'''
row_template = r'     {0:1.1f} & {1:1.1f} & {2:1.1f} & {3:1.1f}  \\'


with open('build/table_Br.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_deg_Br_1, N_Br_1, theta_deg_Br_2, N_Br_2):
        g.write(row_template.format(*row).replace('.',','))
        g.write('\n')
    g.write(table_footer)


# Tabelle für das Absorptionsspektrum für Gallium

theta_deg_Ga_1, theta_deg_Ga_2 = np.array_split(theta_deg_Ga, 2)
N_Ga_1, N_Ga_2 = np.array_split(N_Ga, 2)


table_header = r'''
  \begin{longtable}{S[table-format=2.1] S[table-format=3.1] S[table-format=2.1] S[table-format=3.1]}
  \caption{Messwerte für das Absorptionsspektrum von Gallium}\\
    \toprule
    \multicolumn{1}{c}{ $\theta \:/\:°$ } & \multicolumn{1}{c}{$N\:/\:Imp/s $ }& \multicolumn{1}{c}{ $\theta \:/\:° $}& \multicolumn{1}{c}{ $N\:/\:Imp/s $  }\\
    \cmidrule(lr){1-2} \cmidrule(lr{0,5em}){3-4}
'''
table_footer = r'''    \bottomrule
  \label{tab:4}
  \end{longtable}
'''
row_template = r'     {0:1.1f} & {1:1.1f} & {2:1.1f} & {3:1.1f}  \\'


with open('build/table_Ga.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_deg_Ga_1, N_Ga_1, theta_deg_Ga_2, N_Ga_2):
        g.write(row_template.format(*row).replace('.',','))
        g.write('\n')
    g.write(table_footer)


# Tabelle für das Absorptionsspektrum für Rubidium

theta_deg_Rb_1, theta_deg_Rb_2 = np.array_split(theta_deg_Rb, 2)
N_Rb_1, N_Rb_2 = np.array_split(N_Rb, 2)


table_header = r'''
  \begin{longtable}{S[table-format=2.1] S[table-format=2.1] S[table-format=2.1] S[table-format=2.1]}
  \caption{Messwerte für das Absorptionsspektrum von Rubidium}\\
    \toprule
    \multicolumn{1}{c}{ $\theta \:/\:°$ } & \multicolumn{1}{c}{$N\:/\:Imp/s $ }& \multicolumn{1}{c}{ $\theta \:/\:° $}& \multicolumn{1}{c}{ $N\:/\:Imp/s $  }\\
    \cmidrule(lr){1-2} \cmidrule(lr{0,5em}){3-4}
'''
table_footer = r'''    \bottomrule
  \label{tab:5}
  \end{longtable}
'''
row_template = r'     {0:1.1f} & {1:1.1f} & {2:1.1f} & {3:1.1f}  \\'


with open('build/table_Rb.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_deg_Rb_1, N_Rb_1, theta_deg_Rb_2, N_Rb_2):
        g.write(row_template.format(*row).replace('.',','))
        g.write('\n')
    g.write(table_footer)

# Tabelle für das Absorptionsspektrum für Strontium

theta_deg_Sr_1, theta_deg_Sr_2 = np.array_split(theta_deg_Sr, 2)
N_Sr_1, N_Sr_2 = np.array_split(N_Sr, 2)


table_header = r'''
  \begin{longtable}{S[table-format=2.1] S[table-format=3.1] S[table-format=2.1] S[table-format=3.1]}
  \caption{Messwerte für das Absorptionsspektrum von Strontium}\\
    \toprule
    \multicolumn{1}{c}{ $\theta \:/\:°$ } & \multicolumn{1}{c}{$N\:/\:Imp/s $ }& \multicolumn{1}{c}{ $\theta \:/\:° $}& \multicolumn{1}{c}{ $N\:/\:Imp/s $  }\\
    \cmidrule(lr){1-2} \cmidrule(lr{0,5em}){3-4}
'''
table_footer = r'''    \bottomrule
  \label{tab:6}
  \end{longtable}
'''
row_template = r'     {0:1.1f} & {1:1.1f} & {2:1.1f} & {3:1.1f}  \\'


with open('build/table_Sr.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_deg_Sr_1, N_Sr_1, theta_deg_Sr_2, N_Sr_2):
        g.write(row_template.format(*row).replace('.',','))
        g.write('\n')
    g.write(table_footer)


# Tabelle für das Absorptionsspektrum für Zink

theta_deg_Zn_1, theta_deg_Zn_2 = np.array_split(theta_deg_Zn, 2)
N_Zn_1, N_Zn_2 = np.array_split(N_Zn, 2)


table_header = r'''
  \begin{longtable}{S[table-format=2.1] S[table-format=3.1] S[table-format=2.1] S[table-format=3.1]}
  \caption{Messwerte für das Absorptionsspektrum von Zink}\\
    \toprule
    \multicolumn{1}{c}{ $\theta \:/\:°$ } & \multicolumn{1}{c}{$N\:/\:Imp/s $ }& \multicolumn{1}{c}{ $\theta \:/\:° $}& \multicolumn{1}{c}{ $N\:/\:Imp/s $  }\\
    \cmidrule(lr){1-2} \cmidrule(lr{0,5em}){3-4}
'''
table_footer = r'''    \bottomrule
  \label{tab:7}
  \end{longtable}
'''
row_template = r'     {0:1.1f} & {1:1.1f} & {2:1.1f} & {3:1.1f}  \\'


with open('build/table_Zn.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_deg_Zn_1, N_Zn_1, theta_deg_Zn_2, N_Zn_2):
        g.write(row_template.format(*row).replace('.',','))
        g.write('\n')
    g.write(table_footer)


# Tabelle für das Absorptionsspektrum für Zirkonium

theta_deg_Zr_1, theta_deg_Zr_2 = np.array_split(theta_deg_Zr, 2)
N_Zr_1, N_Zr_2 = np.array_split(N_Zr, 2)


table_header = r'''
  \begin{longtable}{S[table-format=2.1] S[table-format=3.1] S[table-format=2.1] S[table-format=3.1]}
  \caption{Messwerte für das Absorptionsspektrum von Zirkonium}\\
    \toprule
    \multicolumn{1}{c}{ $\theta \:/\:°$ } & \multicolumn{1}{c}{$N\:/\:Imp/s $ }& \multicolumn{1}{c}{ $\theta \:/\:° $}& \multicolumn{1}{c}{ $N\:/\:Imp/s $  }\\
    \cmidrule(lr){1-2} \cmidrule(lr{0,5em}){3-4}
'''
table_footer = r'''    \bottomrule
  \label{tab:8}
  \end{longtable}
'''
row_template = r'     {0:1.1f} & {1:1.1f} & {2:1.1f} & {3:1.1f}  \\'


with open('build/table_Zr.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_deg_Zr_1, N_Zr_1, theta_deg_Zr_2, N_Zr_2):
        g.write(row_template.format(*row).replace('.',','))
        g.write('\n')
    g.write(table_footer)

