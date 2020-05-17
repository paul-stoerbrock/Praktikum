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

theta_bragg, N_bragg = np.genfromtxt('Bragg.dat', unpack = True)
theta_brom, N_brom = np.genfromtxt('Brom.dat', unpack = True)
theta_cu, N_cu = np.genfromtxt('Emissionsspektrum.dat', unpack = True)
theta_gallium, N_gallium = np.genfromtxt('Gallium.dat', unpack = True)
theta_rub, N_rub = np.genfromtxt('Rubidium.dat', unpack = True)
theta_stron, N_stron = np.genfromtxt('Strontium.dat', unpack = True)
theta_zink, N_zink = np.genfromtxt('Zink.dat', unpack = True)
theta_zirk, N_zirk = np.genfromtxt('Zirkonium.dat', unpack = True)

#Globale Variablen---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

d_LiF = 201.4*1e-12 #piko meter
R_infty = 13.6 #eV

#Ordnungszahlen---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

z_cu = 29
z_br = 35
z_ga = 31
z_rb = 37
z_sr = 38
z_zn = 30
z_zr = 40

#Literaturwerte---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

theta_bragg_lit = 28 #degree
E_Kedge_cu = 8987.9615 #eV


#Braggsche Bedingung-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def Gauß(x, a, b, c, d):
  return a*np.exp(-((x-b)/c)**2)+d

popBragg, pcovBragg = curve_fit(
    Gauß,
    theta_bragg[:41],
    N_bragg[:41],
    sigma=None,
    absolute_sigma=True,
    p0=[150, 28, -30, 20]
    )

plt.plot(theta_bragg, N_bragg, 'kx', label='Messwerte')
x_plot = np.linspace(26, 30, 1000)
plt.plot(x_plot, Gauß(x_plot,*popBragg), linestyle='-', label='Nicht-lineare Regression')

plt.xlabel(r'Winkel $\theta\;[\degree]$')
plt.ylabel(r'Impuls [s$^{-1}$]')
ymax_bragg = np.max(N_bragg)
xmax_bragg = theta_bragg[22]
plt.axvline(x=xmax_bragg, color='r', linestyle= '--', label='Maximum')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plotBragg.pdf')
plt.close()

a = popBragg[0]
b = popBragg[1]
c = popBragg[2]
d = popBragg[3]

theta_bragg_diff = np.abs(theta_bragg_lit - xmax_bragg)

# tex file for a

with open('build/a.tex', 'w') as f:
  f.write(make_SI(a,r'', figures=1))

# tex file for b

with open('build/b.tex', 'w') as f:
  f.write(make_SI(b,r'', figures=1))

# tex file for c

with open('build/c.tex', 'w') as f:
  f.write(make_SI(c,r'', figures=1))

# tex file for d

with open('build/d.tex', 'w') as f:
  f.write(make_SI(d,r'', figures=1))

# tex file for theta_bragg_diff

with open('build/theta_bragg_diff.tex', 'w') as f:
  f.write(make_SI(theta_bragg_diff,r'\degree', figures=1))

#Emissionsspektrum Kupfer---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

par, covm = np.polyfit(theta_cu, N_cu, deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(theta_cu, N_cu, 'kx', label='Messwerte')
x_plot = np.linspace(8, 30, 1000)

plt.xlabel(r'Winkel $\theta\;[\degree]$')
plt.ylabel(r'Impuls [s$^{-1}$]')
ymax_cu = np.max(N_cu)
xmax_cu = theta_cu[145]
ylmax_cu = N_cu[122]
xlmax_cu = theta_cu[122]
HWB_lKa_cu_a = 22.3
HWB_rKa_cu_a = 22.9
HWB_lKa_cu_b = 20
HWB_rKa_cu_b = 20.6
plt.axvline(x=xmax_cu, color='r', linestyle= ':', label='$K_{\\alpha}$-Linie')
plt.axvline(x=xlmax_cu, color='b', linestyle= ':', label='$K_{\\beta}$-Linie')
plt.annotate('$K_{\\alpha}$', xy=(xmax_cu, 2500), xycoords='data', xytext=(xmax_cu+1, 2500), textcoords='data', arrowprops=dict(arrowstyle='->', facecolor='grey'), horizontalalignment='left')
plt.annotate('$K_{\\beta}$', xy=(xlmax_cu, 2500), xycoords='data', xytext=(xlmax_cu+1, 2500), textcoords='data', arrowprops=dict(arrowstyle='->', facecolor='grey'), horizontalalignment='left')
plt.annotate('Bremsstrahlung', xy=(xmax_cu/2, 1100), xytext=(xmax_cu/2, 1100), verticalalignment='top', horizontalalignment='left')
plt.axvline(x=HWB_lKa_cu_a, color='grey', linestyle= '-', label='Grenzen der Halbwertsbreite (HWB)')
plt.axvline(x=HWB_rKa_cu_a, color='grey', linestyle= '-', label='')
plt.axvline(x=HWB_lKa_cu_b, color='grey', linestyle= '-', label='')
plt.axvline(x=HWB_rKa_cu_b, color='grey', linestyle= '-', label='')
plt.annotate('$HWB_{K_{\\alpha}}$', xy=(HWB_rKa_cu_a, 2000), xycoords='data', xytext=(HWB_rKa_cu_a+1, 2000), textcoords='data', arrowprops=dict(arrowstyle='->', facecolor='grey'), horizontalalignment='left')
plt.annotate('$HWB_{K_{\\beta}}$', xy=(HWB_lKa_cu_b, 2000), xycoords='data', xytext=(HWB_lKa_cu_b-1, 2000), textcoords='data', arrowprops=dict(arrowstyle='->', facecolor='grey'), horizontalalignment='right')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plotCu.pdf')
plt.close()

xmax_cu = ufloat(np.deg2rad(22.5), np.deg2rad(0.1))
xlmax_cu = ufloat(np.deg2rad(20.2), np.deg2rad(0.1))
Breite_lKa_cu_a = ufloat(np.deg2rad(22.3), np.deg2rad(0.1))
Breite_rKa_cu_a = ufloat(np.deg2rad(22.9), np.deg2rad(0.1))
Breite_lKa_cu_b = ufloat(np.deg2rad(20), np.deg2rad(0.1))
Breite_rKa_cu_b = ufloat(np.deg2rad(20.6), np.deg2rad(0.1))


HWB_cu_bogen_a = abs(Breite_rKa_cu_a-Breite_lKa_cu_a)
HWB_cu_bogen_b = abs(Breite_rKa_cu_b-Breite_lKa_cu_b)


L_cu_a = 2*d_LiF*unp.sin(xmax_cu)
L_cu_b = 2*d_LiF*unp.sin(xlmax_cu)
L_HWB_cu_a = 2*d_LiF*unp.sin(HWB_cu_bogen_a)
L_HWB_cu_b = 2*d_LiF*unp.sin(HWB_cu_bogen_b)


E_cu_a = (const.h*const.c)/(L_cu_a*const.e)
E_cu_b = (const.h*const.c)/(L_cu_b*const.e)
E_HWB_cu_a = (const.h*const.c)/(L_HWB_cu_a*const.e)
E_HWB_cu_b = (const.h*const.c)/(L_HWB_cu_b*const.e)


A_cu_a = E_cu_a/E_HWB_cu_a
A_cu_b = E_cu_b/E_HWB_cu_b


sigma1_cu = z_cu-np.sqrt(E_Kedge_cu/R_infty)
sigma2_cu = z_cu-unp.sqrt(4*(E_Kedge_cu-E_cu_a)/R_infty)
sigma3_cu = z_cu-unp.sqrt(9*(E_Kedge_cu-E_cu_b)/R_infty)

# tex file for xmax_cu

with open('build/xmax_cu.tex', 'w') as f:
  f.write(make_SI(xmax_cu,r'°', figures=1))

# tex file for L_cu_a

with open('build/L_cu_a.tex', 'w') as f:
  f.write(make_SI(L_cu_a,r'', figures=1))

# tex file for E_cu_a

with open('build/E_cu_a.tex', 'w') as f:
  f.write(make_SI(E_cu_a,r'', figures=1))

# tex file for xlmax_cu

with open('build/xlmax_xu.tex', 'w') as f:
  f.write(make_SI(xlmax_cu,r'°', figures=1))

# tex file for L_cu_b

with open('build/L_cu_b.tex', 'w') as f:
  f.write(make_SI(L_cu_b,r'', figures=1))

# tex file for E_cu_b

with open('build/E_cu_b.tex', 'w') as f:
  f.write(make_SI(E_cu_b,r'', figures=1))

# tex file for Breite_lKa_cu_a

with open('build/Breite_lKa_cu_a.tex', 'w') as f:
  f.write(make_SI(Breite_lKa_cu_a,r'°', figures=1))

# tex file for Breite_rKa_cu_a

with open('build/Breite_rKa_cu_a.tex', 'w') as f:
  f.write(make_SI(Breite_rKa_cu_a,r'°', figures=1))

# tex file for Breite_lKa_cu_b

with open('build/Breite_lKa_cu_b.tex', 'w') as f:
  f.write(make_SI(Breite_lKa_cu_b,r'°', figures=1))

# tex file for Breite_rKa_cu_b

with open('build/Breite_rKa_cu_b.tex', 'w') as f:
  f.write(make_SI(Breite_rKa_cu_b,r'°', figures=1))

# tex file for L_HWB_cu_a

with open('build/L_HWB_cu_a.tex', 'w') as f:
  f.write(make_SI(L_HWB_cu_a,r'', figures=1))

# tex file for E_HWB_cu_a

with open('build/E_HWB_cu_a.tex', 'w') as f:
  f.write(make_SI(E_HWB_cu_a,r'', figures=1))

# tex file for A_cu_a

with open('build/A_cu_a.tex', 'w') as f:
  f.write(make_SI(A_cu_a,r'', figures=1))

# tex file for L_HWB_cu_b

with open('build/L_HWB_cu_b.tex', 'w') as f:
  f.write(make_SI(L_HWB_cu_b,r'', figures=1))

# tex file for E_HWB_cu_b

with open('build/E_HWB_cu_b.tex', 'w') as f:
  f.write(make_SI(E_HWB_cu_b,r'', figures=1))

# tex file for A_cu_b

with open('build/A_cu_b.tex', 'w') as f:
  f.write(make_SI(A_cu_b,r'', figures=1))

# tex file for sigma1_cu

with open('build/sigma1_cu.tex', 'w') as f:
  f.write(make_SI(sigma1_cu,r'', figures=1))

# tex file for sigma2_cu

with open('build/sigma2_cu.tex', 'w') as f:
  f.write(make_SI(sigma2_cu,r'', figures=1))

# tex file for sigma3_cu

with open('build/sigma3_cu.tex', 'w') as f:
  f.write(make_SI(sigma3_cu,r'', figures=1))

#Absorbtionsspektrum Brom---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

I_max_br =N_brom[7]
I_min_br =N_brom[2]
I_K_br = I_min_br +(I_max_br-I_min_br)/2
theta_K_br_graph = 13.2
theta_K_br = ufloat(np.deg2rad(13.2), np.deg2rad(0.1))

par, covm = np.polyfit(theta_brom, N_brom, deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(theta_brom, N_brom, 'kx', label='Messwerte')
x_plot = np.linspace(12, 15, 1000)

plt.xlabel(r'Winkel $\theta\;[\degree]$')
plt.ylabel(r'Impuls [s$^{-1}$]')
ymax_br = np.max(N_brom)
xmax_br = theta_brom[7]
xmin_br = theta_brom[2]
plt.axvline(x=xmax_br, color='r', linestyle= '--', label='$K_{Edge}$-Linie (Ende)')
plt.axvline(x=xmin_br, color='b', linestyle= '--', label='$K_{Edge}$-Linie (Anfang)')
plt.axvline(x=theta_K_br_graph, color='grey', linestyle= '--', label='$\\theta$')
plt.axhline(y=I_K_br, color='grey', linestyle= ':', label='$I_K(\\theta)$')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plotBr.pdf')
plt.close()

L_K_br = 2*d_LiF*unp.sin(theta_K_br)

E_K_br = (const.h*const.c)/(L_K_br*const.e)

# tex file for I_max_br

with open('build/I_max_br.tex', 'w') as f:
  f.write(make_SI(I_max_br,r'\per\second', figures=1))

# tex file for I_min_br

with open('build/I_min_br.tex', 'w') as f:
  f.write(make_SI(I_min_br,r'\per\second', figures=1))

# tex file for I_K_br

with open('build/I_K_br.tex', 'w') as f:
  f.write(make_SI(I_K_br,r'\per\second', figures=1))

# tex file for theta_K_br

with open('build/theta_K_br.tex', 'w') as f:
  f.write(make_SI(theta_K_br_graph,r'°', figures=1))

# tex file for L_K_br

with open('build/L_K_br.tex', 'w') as f:
  f.write(make_SI(L_K_br,r'°', figures=1))

# tex file for E_K_br

with open('build/E_K_br.tex', 'w') as f:
  f.write(make_SI(E_K_br,r'°', figures=1))

#Absorbtionsspektrum Gallium-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

I_max_ga =N_gallium[6]
I_min_ga =N_gallium[1]
I_K_ga = I_min_ga +(I_max_ga-I_min_ga)/2
theta_K_ga_graph = 17.3375
theta_K_ga = ufloat(np.deg2rad(17.3375), np.deg2rad(0.1))

par, covm = np.polyfit(theta_gallium, N_gallium, deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(theta_gallium, N_gallium, 'kx', label='Messwerte')
x_plot = np.linspace(17, 20, 1000)

plt.xlabel(r'Winkel $\theta\;[\degree]$')
plt.ylabel(r'Impuls [s$^{-1}$]')
ymax_ga = np.max(N_gallium)
xmax_ga = theta_gallium[6]
xmin_ga = theta_gallium[1]
plt.axvline(x=xmax_ga, color='r', linestyle= '--', label='$K_{Edge}$-Linie (Ende)')
plt.axvline(x=xmin_ga, color='b', linestyle= '--', label='$K_{Edge}$-Linie (Anfang)')
plt.axvline(x=theta_K_ga_graph, color='grey', linestyle= '--', label='$\\theta$')
plt.axhline(y=I_K_ga, color='grey', linestyle= ':', label='$I_K(\\Theta)$')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plotGa.pdf')
plt.close()

L_K_ga = 2*d_LiF*unp.sin(theta_K_ga)

E_K_ga = (const.h*const.c)/(L_K_ga*const.e)

# tex file for I_max_ga

with open('build/I_max_ga.tex', 'w') as f:
  f.write(make_SI(I_max_ga,r'\per\second', figures=1))

# tex file for I_min_ga

with open('build/I_min_ga.tex', 'w') as f:
  f.write(make_SI(I_min_ga,r'\per\second', figures=1))

# tex file for I_K_ga

with open('build/I_K_ga.tex', 'w') as f:
  f.write(make_SI(I_K_ga,r'\per\second', figures=1))

# tex file for theta_K_ga

with open('build/theta_K_ga.tex', 'w') as f:
  f.write(make_SI(theta_K_ga_graph,r'°', figures=1))

# tex file for L_K_ga

with open('build/L_K_ga.tex', 'w') as f:
  f.write(make_SI(L_K_ga,r'°', figures=1))

# tex file for E_K_ga

with open('build/E_K_ga.tex', 'w') as f:
  f.write(make_SI(E_K_ga,r'°', figures=1))

#Absorbtionsspektrum Rubidium-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

I_max_rb =N_rub[9]
I_min_rb =N_rub[2]
I_K_rb = I_min_rb +(I_max_rb-I_min_rb)/2
theta_K_rb_graph = 11.774
theta_K_rb = ufloat(np.deg2rad(11.774), np.deg2rad(0.1))

par, covm = np.polyfit(theta_rub, N_rub, deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(theta_rub, N_rub, 'kx', label='Messwerte')
x_plot = np.linspace(11, 13, 1000)

plt.xlabel(r'Winkel $\theta\;[\degree]$')
plt.ylabel(r'Impuls [s$^{-1}$]')
ymax_rb = np.max(N_rub)
xmax_rb = theta_rub[9]
xmin_rb = theta_rub[2]
plt.axvline(x=xmax_rb, color='r', linestyle= '--', label='$K_{Edge}$-Linie (Ende)')
plt.axvline(x=xmin_rb, color='b', linestyle= '--', label='$K_{Edge}$-Linie (Anfang)')
plt.axvline(x=theta_K_rb_graph, color='grey', linestyle= '--', label='$\\theta$')
plt.axhline(y=I_K_rb, color='grey', linestyle= ':', label='$I_K(\\Theta)$')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plotRb.pdf')
plt.close()

L_K_rb = 2*d_LiF*unp.sin(theta_K_rb)

E_K_rb = (const.h*const.c)/(L_K_rb*const.e)

# tex file for I_max_rb

with open('build/I_max_rb.tex', 'w') as f:
  f.write(make_SI(I_max_rb,r'\per\second', figures=1))

# tex file for I_min_rb

with open('build/I_min_rb.tex', 'w') as f:
  f.write(make_SI(I_min_rb,r'\per\second', figures=1))

# tex file for I_K_rb

with open('build/I_K_rb.tex', 'w') as f:
  f.write(make_SI(I_K_rb,r'\per\second', figures=1))

# tex file for theta_K_rb

with open('build/theta_K_br.tex', 'w') as f:
  f.write(make_SI(theta_K_rb_graph,r'°', figures=1))

# tex file for L_K_rb

with open('build/L_K_rb.tex', 'w') as f:
  f.write(make_SI(L_K_rb,r'°', figures=1))

# tex file for E_K_rb

with open('build/E_K_rb.tex', 'w') as f:
  f.write(make_SI(E_K_rb,r'°', figures=1))

#Absorbtionsspektrum Strontium-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

I_max_sr =N_stron[9]
I_min_sr =N_stron[2]
I_K_sr = I_min_sr +(I_max_sr-I_min_sr)/2
theta_K_sr_graph = 11.087
theta_K_sr = ufloat(np.deg2rad(11.087), np.deg2rad(0.1))

par, covm = np.polyfit(theta_stron, N_stron, deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(theta_stron, N_stron, 'kx', label='Messwerte')
x_plot = np.linspace(10, 13, 1000)

plt.xlabel(r'Winkel $\theta\;[\degree]$')
plt.ylabel(r'Impuls [s$^{-1}$]')
ymax_sr = np.max(N_stron)
xmax_sr = theta_stron[9]
xmin_sr = theta_stron[2]
plt.axvline(x=xmax_sr, color='r', linestyle= '--', label='$K_{Edge}$-Linie (Ende)')
plt.axvline(x=xmin_sr, color='b', linestyle= '--', label='$K_{Edge}$-Linie (Anfang)')
plt.axvline(x=theta_K_sr_graph, color='grey', linestyle= '--', label='$\\theta$')
plt.axhline(y=I_K_sr, color='grey', linestyle= ':', label='$I_K(\\Theta)$')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plotSr.pdf')
plt.close()

L_K_sr = 2*d_LiF*unp.sin(theta_K_sr)

E_K_sr = (const.h*const.c)/(L_K_sr*const.e)

# tex file for I_max_sr

with open('build/I_max_sr.tex', 'w') as f:
  f.write(make_SI(I_max_sr,r'\per\second', figures=1))

# tex file for I_min_sr

with open('build/I_min_sr.tex', 'w') as f:
  f.write(make_SI(I_min_sr,r'\per\second', figures=1))

# tex file for I_K_sr

with open('build/I_K_sr.tex', 'w') as f:
  f.write(make_SI(I_K_sr,r'\per\second', figures=1))

# tex file for theta_K_sr

with open('build/theta_K_sr.tex', 'w') as f:
  f.write(make_SI(theta_K_sr_graph,r'°', figures=1))

# tex file for L_K_sr

with open('build/L_K_sr.tex', 'w') as f:
  f.write(make_SI(L_K_sr,r'°', figures=1))

# tex file for E_K_sr

with open('build/E_K_sr.tex', 'w') as f:
  f.write(make_SI(E_K_sr,r'°', figures=1))

#Absorbtionsspektrum Zink-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

I_max_zn =N_zink[10]
I_min_zn =N_zink[4]
I_K_zn = I_min_zn +(I_max_zn-I_min_zn)/2
theta_K_zn_graph = 18.67
theta_K_zn = ufloat(np.deg2rad(18.67), np.deg2rad(0.1))

par, covm = np.polyfit(theta_zink, N_zink, deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(theta_zink, N_zink, 'kx', label='Messwerte')
x_plot = np.linspace(18, 20, 1000)

plt.xlabel(r'Winkel $\theta\;[\degree]$')
plt.ylabel(r'Impuls [s$^{-1}$]')
ymax_zn = np.max(N_zink)
xmax_zn = theta_zink[10]
xmin_zn = theta_zink[4]
plt.axvline(x=xmax_zn, color='r', linestyle= '--', label='$K_{Edge}$-Linie (Ende)')
plt.axvline(x=xmin_zn, color='b', linestyle= '--', label='$K_{Edge}$-Linie (Anfang)')
plt.axvline(x=theta_K_zn_graph, color='grey', linestyle= '--', label='$\\theta$')
plt.axhline(y=I_K_zn, color='grey', linestyle= ':', label='$I_K(\\Theta)$')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plotZn.pdf')
plt.close()

L_K_zn = 2*d_LiF*unp.sin(theta_K_zn)

E_K_zn = (const.h*const.c)/(L_K_zn*const.e)

# tex file for I_max_zn

with open('build/I_max_zn.tex', 'w') as f:
  f.write(make_SI(I_max_zn,r'\per\second', figures=1))

# tex file for I_min_zn

with open('build/I_min_zn.tex', 'w') as f:
  f.write(make_SI(I_min_zn,r'\per\second', figures=1))

# tex file for I_K_zn

with open('build/I_K_zn.tex', 'w') as f:
  f.write(make_SI(I_K_zn,r'\per\second', figures=1))

# tex file for theta_K_zn

with open('build/theta_K_zn.tex', 'w') as f:
  f.write(make_SI(theta_K_zn_graph,r'°', figures=1))

# tex file for L_K_zn

with open('build/L_K_zn.tex', 'w') as f:
  f.write(make_SI(L_K_zn,r'°', figures=1))

# tex file for E_K_zn

with open('build/E_K_zn.tex', 'w') as f:
  f.write(make_SI(E_K_zn,r'°', figures=1))

#Absorbtionsspektrum Zirkonium-------------------------------------------------------------------------------------------------------------------------------------------------------------------------

I_max_zr =N_zirk[9]
I_min_zr =N_zirk[0]
I_K_zr = I_min_zr +(I_max_zr-I_min_zr)/2
theta_K_zr_graph = 9.96
theta_K_zr = ufloat(np.deg2rad(9.96), np.deg2rad(0.1))

par, covm = np.polyfit(theta_zirk, N_zirk, deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(theta_zirk, N_zirk, 'kx', label='Messwerte')
x_plot = np.linspace(9, 11, 1000)

plt.xlabel(r'Winkel $\theta\;[\degree]$')
plt.ylabel(r'Impuls [s$^{-1}$]')
ymax_zr = np.max(N_zirk)
xmax_zr = theta_zirk[9]
xmin_zr = theta_zirk[0]
plt.axvline(x=xmax_zr, color='r', linestyle= '--', label='$K_{Edge}$-Linie (Ende)')
plt.axvline(x=xmin_zr, color='b', linestyle= '--', label='$K_{Edge}$-Linie (Anfang)')
plt.axvline(x=theta_K_zr_graph, color='grey', linestyle= '--', label='$\\theta$')
plt.axhline(y=I_K_zr, color='grey', linestyle= ':', label='$I_K(\\Theta)$')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plotZr.pdf')
plt.close()

L_K_zr = 2*d_LiF*unp.sin(theta_K_zr)

E_K_zr = (const.h*const.c)/(L_K_zr*const.e)

# tex file for I_max_zr

with open('build/I_max_zr.tex', 'w') as f:
  f.write(make_SI(I_max_zr,r'\per\second', figures=1))

# tex file for I_min_zr

with open('build/I_min_zr.tex', 'w') as f:
  f.write(make_SI(I_min_zr,r'\per\second', figures=1))

# tex file for I_K_zr

with open('build/I_K_zr.tex', 'w') as f:
  f.write(make_SI(I_K_zr,r'\per\second', figures=1))

# tex file for theta_K_zr

with open('build/theta_K_zr.tex', 'w') as f:
  f.write(make_SI(theta_K_zr_graph,r'°', figures=1))

# tex file for L_K_zr

with open('build/L_K_zr.tex', 'w') as f:
  f.write(make_SI(L_K_zr,r'°', figures=1))

# tex file for E_K_zr

with open('build/E_K_zr.tex', 'w') as f:
  f.write(make_SI(E_K_zr,r'°', figures=1))

#Rydbergfrequenz R_infty-------------------------------------------------------------------------------------------------------------------------------------------------

sigma_br = z_br - unp.sqrt(E_K_br/R_infty-const.alpha**2*z_br**4/4)
sigma_ga = z_ga - unp.sqrt(E_K_ga/R_infty-const.alpha**2*z_ga**4/4)
sigma_rb = z_rb - unp.sqrt(E_K_rb/R_infty-const.alpha**2*z_rb**4/4)
sigma_sr = z_sr - unp.sqrt(E_K_sr/R_infty-const.alpha**2*z_sr**4/4)
sigma_zn = z_zn - unp.sqrt(E_K_zn/R_infty-const.alpha**2*z_zn**4/4)
sigma_zr = z_zr - unp.sqrt(E_K_zr/R_infty-const.alpha**2*z_zr**4/4)

# tex file for sigma_br

with open('build/sigma_br.tex', 'w') as f:
  f.write(make_SI(sigma_br,r'', figures=1))

# tex file for sigma_ga

with open('build/sigma_ga.tex', 'w') as f:
  f.write(make_SI(sigma_ga,r'', figures=1))

# tex file for sigma_rb

with open('build/sigma_rb.tex', 'w') as f:
  f.write(make_SI(sigma_rb,r'', figures=1))

# tex file for sigma_sr

with open('build/sigma_sr.tex', 'w') as f:
  f.write(make_SI(sigma_sr,r'', figures=1))

# tex file for sigma_zn

with open('build/sigma_zn.tex', 'w') as f:
  f.write(make_SI(sigma_zn,r'', figures=1))

# tex file for sigma_zr

with open('build/sigma_zr.tex', 'w') as f:
  f.write(make_SI(sigma_zr,r'', figures=1))

sigma = np.array([sigma_br.n, sigma_ga.n, sigma_rb.n, sigma_sr.n, sigma_zn.n, sigma_zr.n])
z = np.array([z_br, z_ga, z_rb, z_sr, z_zn, z_zr])
z_eff = z - sigma

E_K = np.array([E_K_br.n, E_K_ga.n, E_K_rb.n, E_K_sr.n, E_K_zn.n, E_K_zr.n])
E_K_sqrt = np.sqrt(E_K)

par, covm = np.polyfit(z_eff, E_K_sqrt, deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(z_eff, E_K_sqrt, 'kx', label='Messwerte')
x_plot = np.linspace(26, 36, 1000)
plt.plot(x_plot ,par[0]*x_plot+par[1], 'b-', label="Lineare Regression")
plt.legend(loc="best")
plt.xlabel(r'$Z_{eff}$')
plt.ylabel(r'$\sqrt{E_K}$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotE_K.pdf')
plt.close()

m_ry = ufloat(par[0], err[0])
b_ry = ufloat(par[1], err[1])

# tex file for m_ry

with open('build/m_ry.tex', 'w') as f:
  f.write(make_SI(m_ry,r'', figures=1))

# tex file for b_ry

with open('build/b_ry.tex', 'w') as f:
  f.write(make_SI(b_ry,r'', figures=1))

Rydb = m_ry**2
print(Rydb)
# tex file for Rydb

with open('build/Rydb.tex', 'w') as f:
  f.write(make_SI(Rydb,r'', figures=1))

#Tabelle Kupfer-------------------------------------------------------------------------------------------------------------------------------------------------

theta_cu1, theta_cu2, theta_cu3 = np.array_split(theta_cu,3)
N_Cu1, N_Cu2, N_Cu3 = np.array_split(N_cu,3)

table_header = r'''
\begin{longtable}{S[table-format=2.1] S[table-format=3.0] S[table-format=4.1] S[table-format=2.0] S[table-format=3.1] S[table-format=4.0]}
    \caption{Impulsrate $N$ der Kupferröhre in Abh\"angigkeit des Streuwinkels $\theta$}\\
    \toprule
    {$\theta \; [\si{\degree}]$} & {$\text{Impulse} \; N_{\text{Kupfer}} \; [\si{\per\second}]$} & 
    {$\theta \; [\si{\degree}]$} & {$\text{Impulse} \; N_{\text{Kupfer}} \; [\si{\per\second}]$} & 
    {$\theta \; [\si{\degree}]$} & {$\text{Impulse} \; N_{\text{Kupfer}} \; [\si{\per\second}]$}\\
    \cmidrule(lr{0.5em}){1-2} \cmidrule(lr{0.5em}){3-4} \cmidrule(lr{0.5em}){5-6}
'''
table_footer = r''' 
    \bottomrule
    \label{tab:1}
\end{longtable}
'''

row_template = r'     {0:2.1f} & {1:3.0f} & {2:4.1f} & {3:2.0f} & {4:3.1f} & {5:4.0f}\\'

with open('build/table_Cu.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_cu1, N_Cu1, theta_cu2, N_Cu2, theta_cu3, N_Cu3):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)


#Tabelle Bragg-------------------------------------------------------------------------------------------------------------------------------------------------

theta_bragg1, theta_bragg2 = np.array_split(theta_bragg, 2)
N_bragg1, N_bragg2 = np.array_split(N_bragg, 2)

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    {$\text{Winkel} \; \theta\; [\si{\degree}]$} & {$\text{Impulse $N_{\text{Bragg}}$} \; [\si{\per\second}]$} &
    {$\text{Winkel} \; \theta\; [\si{\degree}]$} & {$\text{Impulse $N_{\text{Bragg}}$} \; [\si{\per\second}]$} \\

    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f} & {2:1.2f} & {3:1.2f}  \\'

with open('build/table_bragg.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_bragg1, N_bragg1, theta_bragg2, N_bragg2):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

#Tabelle Brom-------------------------------------------------------------------------------------------------------------------------------------------------

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    {$\text{Winkel} \; \theta\; [\si{\degree}]$} & {$\text{Impulse $N_{\text{Brom}}$} \; [\si{\per\second}]$} \\

    \cmidrule(lr{0,5em}){1-2}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f}  \\'

with open('build/table_brom.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_brom, N_brom):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

#Tabelle Gallium-------------------------------------------------------------------------------------------------------------------------------------------------

theta_gallium1, theta_gallium2 = np.array_split(theta_gallium, 2)
N_gallium1, N_gallium2 = np.array_split(N_gallium, 2)

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    {$\text{Winkel} \; \theta\; [\si{\degree}]$} & {$\text{Impulse $N_{\text{Gallium}}$} \; [\si{\per\second}]$} &
    {$\text{Winkel} \; \theta\; [\si{\degree}]$} & {$\text{Impulse $N_{\text{Gallium}}$} \; [\si{\per\second}]$} \\

    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f} & {2:1.2f} & {3:1.2f}  \\'

with open('build/table_gallium.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_gallium1, N_gallium1, theta_gallium2, N_gallium2):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

#Tabelle Rubidium-------------------------------------------------------------------------------------------------------------------------------------------------

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    {$\text{Winkel} \; \theta\; [\si{\degree}]$} & {$\text{Impulse $N_{\text{Rubidium}}$} \; [\si{\per\second}]$} \\

    \cmidrule(lr{0,5em}){1-2}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f}  \\'

with open('build/table_rub.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_rub, N_rub):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

#Tabelle Rubidium-------------------------------------------------------------------------------------------------------------------------------------------------

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    {$\text{Winkel} \; \theta\; [\si{\degree}]$} & {$\text{Impulse $N_{\text{Strontium}}$} \; [\si{\per\second}]$} \\

    \cmidrule(lr{0,5em}){1-2}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f}  \\'

with open('build/table_stron.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_stron, N_stron):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

#Tabelle Zink-------------------------------------------------------------------------------------------------------------------------------------------------

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    {$\text{Winkel} \; \theta\; [\si{\degree}]$} & {$\text{Impulse $N_{\text{Zink}}$} \; [\si{\per\second}]$} \\

    \cmidrule(lr{0,5em}){1-2}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f}  \\'

with open('build/table_zink.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_zink, N_zink):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

#Tabelle Zirkonium-------------------------------------------------------------------------------------------------------------------------------------------------

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    {$\text{Winkel} \; \theta\; [\si{\degree}]$} & {$\text{Impulse $N_{\text{Zirkonium}}$} \; [\si{\per\second}]$} \\

    \cmidrule(lr{0,5em}){1-2}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f}  \\'

with open('build/table_zirk.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_zirk, N_zirk):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)