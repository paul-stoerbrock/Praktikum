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


# Funktionsdefinitionen #########################################################################



# Messdaten ##################################################################################
Aa, ta = np.genfromtxt('4a.txt', unpack=True) # Aa=Amplitudenspannung für 4a), ta=Zeitdifferenz der Nulldurchgänge für 4a)

f, A, t = np.genfromtxt('data.txt', unpack=True) # f=Frequenz, A=Amplitudenspannung, t=Zeitdifferenz der Nulldurhgänge

U0 = 10 # angelegte Spannung in Volt

Umax= 13  # Resonanzamplitude

A0 = A/U0

phi = f * t * 2 * np.pi # Umrechnung von t in phi

R= 2.9 * 1e03      # Widerstand für aperiodischen Grenzfall

nuGraph = 14.9*1e03 # nu+ - nu- Breite der Resonanzkurve

nu1 = 3.2 

nu2 = 4.43

nures = 3.74

# Literaturwerte ##########################################################################################

L= ufloat(3.5,0.01)*1e-03 
R1= ufloat(30.3, 0.1)
R2= ufloat(271.6,0.2)
C= 5.01 * 1e-09

# Berechnung R_ap aus Literaturwerten

R_aplit= unp.sqrt(4*L/C)

# Relativer Fehler R

RF_R = (R-R_aplit)/R_aplit

# Berechnung von nu+ - nu-

nuRech= (R2)/(2 * np.pi * L)

# Relativer Fehler nu+ - nu-

RF_nu= (nuGraph-nuRech.n)/nuRech.n



#Plot für a) #################################################################################

slope, intercept, r_value, p_value, std_err = stats.linregress(ta, np.log(Aa))

plt.plot(ta*1e+04,np.log(Aa), 'kx', label="Messdaten")
x_plot = np.linspace(0, 2e-04, 10)
plt.plot(x_plot*1e04,intercept+slope*x_plot, 'r-', label="Lineare Regression")
plt.legend(loc="best")
plt.xlabel('t / s')
plt.ylabel(r'$\log (A) \:/\: $V')
plt.grid()
plt.tight_layout
plt.savefig('build/plota.pdf')
plt.close()

#Plot für c) #######################################################################################################

# Erstellung des Plots c) in ln-Darstellung

plt.plot(f, A0, 'rx', label="Messdaten")
plt.xscale('log')
plt.xlabel(r'$\nu \:/\: $Hz')
plt.ylabel(r'$U_C/U_0 \:/\: $V')
plt.grid()
plt.tight_layout
plt.savefig('build/plotcln.pdf')
plt.close()

# Erstellung des Plots c) in linearer Darstellung

plt.plot(f*1e-03, A0, 'rx', label="Messdaten")
plt.legend(loc="best")
plt.axhline(Umax/(U0*np.sqrt(2)), linestyle='--') 
plt.axvline(27.8,linestyle='--') # nu-
plt.axvline(42.7,linestyle='--') # nu+
plt.xlabel(r'$\nu \:/\: $kHz')
plt.ylabel(r'$U_C/U_0 \:/\: $V')
plt.grid()
plt.tight_layout
plt.savefig('build/plotclin.pdf')
plt.close()


# Plot für d) ###########################################################################################

# Erstellung des Plots d) in ln-Darstellung

plt.plot(f*1e-04, phi, 'kx', label="Messdaten")
plt.xscale('log')
plt.yticks( [0,np.pi/4 ,np.pi/2,3 * np.pi/4 ,np.pi],
            [r'$0$' ,r'$+\pi/4$' ,r'$+\pi/2$' ,r'$+3\pi/4$'  ,r'$+\pi$']
    )
plt.legend(loc="best")
plt.xlabel(r'$\nu \:/\: 10^4$Hz')
plt.ylabel(r'$\phi \:/\: $rad')
plt.tight_layout
plt.grid()
plt.savefig('build/plotdln.pdf')
plt.close()

#Erstellung des Plots d) in linearer Darstellung

plt.plot(f*1e-04, phi, 'kx', label="Messdaten")
plt.xscale('log')
plt.yticks( [0,np.pi/4 ,np.pi/2,3 * np.pi/4 ,np.pi],
            [r'$0$' ,r'$+\pi/4$' ,r'$+\pi/2$' ,r'$+3\pi/4$'  ,r'$+\pi$']
    )
plt.xticks( [3 ,4 ],
            [3,4]
)
plt.legend(loc="best")
plt.xlabel(r'$\nu \:/\: 10^4$ Hz')
plt.ylabel(r'$\phi \:/\: $rad')
plt.axvline(3.2,linestyle='--')   # nu1
plt.axvline(3.74,linestyle='--')  # nures
plt.axvline(4.43,linestyle='--')  # nu2
plt.grid()
plt.tight_layout
plt.savefig('build/plotdlin.pdf')
plt.close()


# SI-Einheiten ####################################################################################






# Erstellung Tabelle c), d) ###################################################################################

f1, f2, f3= np.array_split(f*1e-03, 3)
A1, A2, A3 = np.array_split(A, 3)
t1, t2, t3 = np.array_split(t, 3)


table_header = r'''
  \begin{tabular}{c c c c c c c c c}
    \toprule
    {$\nu \:/\: \si{\kilo\hertz}$}& {$U_C \:/\: \si{\volt}$} & {$t \:/\: \si{\micro\second}$} &
    {$\nu \:/\: \si{\kilo\hertz}$}& {$U_C \:/\: \si{\volt}$} & {$t \:/\: \si{\micro\second}$} &
    {$\nu \:/\: \si{\kilo\hertz}$}& {$U_C \:/\: \si{\volt}$} & {$t \:/\: \si{\micro\second}$} \\


    \cmidrule(lr{0,5em}){1-3} \cmidrule(lr{0,5em}){4-6} \cmidrule(lr{0,5em}){7-9}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f} & {2:1.2f} & {0:1.2f} & {1:1.2f} & {2:1.2f} & {0:1.2f} & {1:1.2f} & {2:1.2f}  \\'

# Tabelle für 4a) wird im Tex Format geschrieben ############################################################################################################################################################################

with open('build/table_c.tex', 'w') as h:
    h.write(table_header)
    for row in zip(f1, A1, t1, f2, A2, t2, f3, A3, t3):
        h.write(row_template.format(*row))
        h.write('\n')
    h.write(table_footer)



# Kontrollprints

print(R_aplit)
print(RF_R)
print(np.exp(1)**intercept)
print(nuRech)
print(nuGraph)
print(RF_nu)