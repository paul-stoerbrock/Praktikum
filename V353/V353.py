
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

# A(Omega)/U0 ############################################################################################################################################################################

def g(f, RC):
    return  1/(np.sqrt(4 * (np.pi)**2 * f**2 * RC**2+1))

def h(f, RC):
    return np.arctan(-2 * np.pi * f * RC)

def d(A0):
    return np.arccos(A0)



U0 = 0.3 

U, tc = np.genfromtxt('entladung.txt', unpack=True) #Variablen definieren U=Q/C, tc=Zeit aus 4a)
U *= 1/5
Uc = np.log(U/U0)
f, A, t = np.genfromtxt('data.txt', unpack=True) #Variablen definieren f=Frequenz, A=Amplitude, t=Zeitabstand der Nulldurchgänge
A *= 1/5
A0 = A / U0
phi = f * t * 2 * np.pi

# Bestimmung von RC (slope, popt, phiRC) ############################################################################################################################################################################

slope, intercept, r_value, p_value, std_err = stats.linregress(tc, np.log(U/U0))

x = ufloat(slope, std_err)
y = -1/x

R = ufloat(15.056e+03, 0.6e+03)+600
C = 93.2e-09
L = R*C



popt, pcov = curve_fit(
    g,
    f,
    A0*5,
    sigma=None,
    absolute_sigma=True,
    p0=[1e-03]
    )

phiRC, phicov = curve_fit(
    h,
    f,
    phi,
    sigma=None,
    absolute_sigma=True,
    p0=[1e-03]
    )

SF2 = (popt[0]-L.n)/L.n
SF1 = (y.n * 1e-03 - L.n)/L.n
SF3 = (phiRC[0] * 1e-03 - L.n)/L.n

# Plot für 4a-d) [A(w)/U0] ############################################################################################################################################################################

plt.plot(tc, U/U0, 'kx', label="Messdaten")
plt.yscale('log')
plt.plot(tc, np.exp(intercept + slope*tc), 'r-', label="Lineare Regression")
plt.legend(loc="best")
plt.xlabel('Zeit in $ms$')
plt.ylabel('$U_c/U_0$ in Volt')
plt.tight_layout
plt.savefig('build/plot4a.pdf')
plt.close()

plt.plot(f, A0*5, 'kx', label="Messdaten")
plt.xscale('log')
x_plot = np.linspace(1, 100000, 100000)
plt.plot(x_plot, g(x_plot,*popt), 'r-', label="Nichtlineare Regression")
plt.legend(loc="best")
plt.xlabel('Frequenz in Hertz')
plt.ylabel('A/$U_0$ in Volt')
plt.grid(True)
plt.tight_layout
plt.savefig('build/plot4b.pdf')
plt.close()

plt.plot(f, phi, 'kx', label="Messdaten")
plt.xscale('log')
x_plot = np.linspace(1, 100000, 100000)
plt.plot(x_plot, h(x_plot,*phiRC), 'r-', label="Nichtlineare Regression")
plt.legend(loc="best")
plt.xlabel('Frequenz in Hertz')
plt.ylabel('Phase in Bogenmaß')
plt.grid(True)
plt.tight_layout
plt.savefig('build/plot4c.pdf')
plt.close()

plt.polar(phi, A0*5, 'kx', label="Messdaten")
x_plot = np.linspace(0, np.pi, 10000)
plt.polar(d(x_plot), x_plot, 'r-', label="Polarplot")
plt.legend(loc="best")
plt.grid(True)
plt.tight_layout
plt.savefig('build/plot4d.pdf')
plt.close()


#richtiger Plot 4c) #########################################################################################################

fr, Ar, tr = np.genfromtxt('true.txt', unpack=True)
Ar *= 1/5

phir = fr * tr * 2 * np.pi

phiRCr, phicovr = curve_fit(
    h,
    fr,
    phir,
    sigma=None,
    absolute_sigma=True,
    p0=[1e-03]
    )

plt.plot(fr, phir, 'kx', label="Messdaten")
plt.xscale('log')
x_plot = np.linspace(1, 100000, 100000)
plt.plot(x_plot, h(x_plot,*phiRCr), 'r-', label="Nichtlineare Regression")
plt.legend(loc="best")
plt.xlabel('Frequenz in Hertz')
plt.ylabel('Phase in Bogenmaß')
plt.grid(True)
plt.tight_layout
plt.savefig('build/plot4ctrue.pdf')
plt.close()

######################################################################################################################################

# richtiger Plot für b) ######################################################################################################

popttrue, pcovtrue = curve_fit(
    g,
    f,
    A0,
    sigma=None,
    absolute_sigma=True,
    p0=[1e-03]
    )

plt.plot(f, A0, 'kx', label="Messdaten")
plt.xscale('log')
x_plot = np.linspace(1, 100000, 100000)
plt.plot(x_plot, g(x_plot,*popttrue), 'r-', label="Nichtlineare Regression")
plt.legend(loc="best")
plt.xlabel('Frequenz in Hertz')
plt.ylabel('A/$U_0$ in Volt')
plt.grid(True)
plt.tight_layout
plt.savefig('build/plot4btrue.pdf')
plt.close()

# richtiger Plot für d) ######################################################################################################

plt.polar(phir, A0, 'kx', label="Messdaten")
x_plot = np.linspace(0, 1, 10000)
plt.polar(d(x_plot), x_plot, 'r-', label="Polarplot")
plt.legend(loc="best")
plt.grid(True)
plt.tight_layout
plt.savefig('build/plot4dtrue.pdf')
plt.close()

#SI Einheiten 

ohmF=r' }{\ohm\farad}$'

#slope und intercept übergeben ##################################################################################################################

with open('build/slope.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{slope:.2f}\pm {std_err:.2f} e-3') 
    RC.write(ohmF)

with open('build/intercept.tex', 'w') as RC:
    RC.write('$\\num{')
    RC.write(f'{intercept:.2f}')
    RC.write('}$')

#RC nach a berechnet ####################################################################################################################

with open('build/mean_aRC.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{y.n:.2f} \pm {y.s:.2f} e-3')
    RC.write(ohmF)

#RC nach b berechnet ########################################################################################################

with open('build/mean_bRC.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{(popt[0]*1e+03):.2f}e-03')
    RC.write(ohmF)

#Literaturwert L nach b berechnet ########################################################################################################

with open('build/L.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{L.n*1000:.3f} \pm {L.s*1000:.3f} e-3')
    RC.write(ohmF)

#Relativen Fehler nach b berechnet ########################################################################################################

with open('build/SF1.tex', 'w') as RC:
    RC.write('$\\num{')
    RC.write(f'{SF1:.3f}')
    RC.write('}$')

with open('build/SF2.tex', 'w') as RC:
    RC.write('$\\num{')
    RC.write(f'{SF2:.3f}')
    RC.write('}$')

#RC nach c berechnet ###################################################################################

with open('build/mean_cRC.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{phiRC[0]:.2f}')
    RC.write(ohmF)

#Relativen Fehler nach b berechnet ########################################################################################################

with open('build/SF3.tex', 'w') as RC:
    RC.write('$\\num{')
    RC.write(f'{SF3:.3f}')
    RC.write('}$')

# Tabelle für 4a) wird erstellt ############################################################################################################################################################################

U1, U2, U3= np.array_split(U, 3)
tc1, tc2, tc3 = np.array_split(tc, 3)

table_header = r'''
  \begin{tabular}{c c c c c c}
    \toprule
    {$U_C \:/\: \si{\volt}$} & {$t \:/\: \si{\milli\second}$} & 
    {$U_C \:/\: \si{\volt}$} & {$t \:/\: \si{\milli\second}$} &  
    {$U_C \:/\: \si{\volt}$} & {$t \:/\: \si{\milli\second}$}\\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4} \cmidrule(lr{0,5em}){5-6}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f} & {2:1.2f} & {3:1.2f} & {4:1.2f} & {5:1.2f}\\'

# Tabelle für 4a) wird im Tex Format geschrieben ############################################################################################################################################################################

with open('build/table_4a.tex', 'w') as h:
    h.write(table_header)
    for row in zip(U1, tc1, U2, tc2, U3, tc3):
        h.write(row_template.format(*row))
        h.write('\n')
    h.write(table_footer)

# Tabelle für 4b-c) wird erstellt ############################################################################################################################################################################

f1, f2, f3= np.array_split(f, 3)
A1, A2, A3 = np.array_split(A, 3)
t1, t2, t3 = np.array_split(1000*t, 3)

table_header = r'''
  \begin{tabular}{c c c c c c c c c}
    \toprule
    {$\nu \:/\: \si{\hertz}$} & {$A(\omega) \:/\: \si{\milli\volt}$} & {$\Delta T \:/\: \si{\milli\second}$} & 
    {$\nu \:/\: \si{\hertz}$} & {$A(\omega) \:/\: \si{\milli\volt}$} & {$\Delta T \:/\: \si{\milli\second}$} & 
    {$\nu \:/\: \si{\hertz}$} & {$A(\omega) \:/\: \si{\milli\volt}$} & {$\Delta T \:/\: \si{\micro\second}$}\\
    \cmidrule(lr{0,5em}){1-3} \cmidrule(lr{0,5em}){4-6} \cmidrule(lr{0,5em}){7-9}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.0f} & {1:1.2f} & {2:1.4f} & {3:1.0f} & {4:1.2f} & {5:1.4f} & {6:1.0f} & {7:1.2f} & {8:1.4f} \\'

# Tabelle für 4b-c) wird im Tex Format geschrieben ############################################################################################################################################################################

with open('build/table_4b.tex', 'w') as i:
    i.write(table_header)
    for row in zip(f1, A1*1e+03, t1, f2, A2*1e+03, t2, f3, A3*1e+03, t3):
        i.write(row_template.format(*row))
        i.write('\n')
    i.write(table_footer)


# Testprints #########################################################################

print(A0)
print(L.n)
print(y.n-L.n)
print(SF1)
print(SF2)