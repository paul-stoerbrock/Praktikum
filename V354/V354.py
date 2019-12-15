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

#Funktionen


#Auswertungsfunktionen

def U(f, R, L, C):
 return 1/(np.sqrt((1-L*C*4*np.pi**2*f**2)**2+4*np.pi**2*f**2*R**2*C**2))

def phi1(f, R, L, C):
    return np.arctan((-2*np.pi*f*R*C)/(1-L*C*4 *np.pi**2*f**2))

# Messdaten ##################################################################################

Aa, ta = np.genfromtxt('4a.txt', unpack=True) # Aa=Amplitudenspannung für 4a), ta=Zeitdifferenz der Nulldurchgänge für 4a)

f, A, t = np.genfromtxt('data.txt', unpack=True) # f=Frequenz, A=Amplitudenspannung, t=Zeitdifferenz der Nulldurhgänge

U0 = 10 # angelegte Spannung in Volt

Umax= 13  # Resonanzamplitude in Volt

A0 = A/U0

phi = f * t * 2 * np.pi # Umrechnung von t in phi

R= 2.9 * 1e03      # Widerstand für aperiodischen Grenzfall

nuGraph = 14.9 *1e03 *2*np.pi # nu+ - nu- Breite der Resonanzkurve

nu1 = 3.2 * 1e04

nu2 = 4.43 * 1e04

nu_res = 3.74 * 1e04

# Literaturwerte ##########################################################################################

L= ufloat(3.5,0.01)*1e-03 
R1= ufloat(30.3, 0.1)
R2= ufloat(271.6,0.2)
C= 5.01 *1e-09

# Berechnung R_ap aus Literaturwerten

R_ap_lit = unp.sqrt(4*L/C)/2*2

# Relativer Fehler R

RF_R = (R-R_ap_lit.n)/R_ap_lit.n # RF= relativer Fehler

# Berechnung von nu+ - nu-

nuRech= (R2+2.5)/(L)

# Relativer Fehler nu+ - nu-

RF_nu= (nuGraph-nuRech.n)/nuRech.n

# Berechnung Omega_0

w0 = 35.3*1e+03*2*np.pi

# Berechnung von w0_Lit.

w0_lit = unp.sqrt(1/(L*C))

# Berechnung von q -> w0/(nu+ - nu-)

q = w0/nuGraph

# Berechnung von q_lit

q_lit = w0_lit/nuRech

# Berechnung des relativen Fehlers von q

RF_q = (q-q_lit.n)/q_lit.n

# Berechnung nu 1 starke Dämfung

nu1_lit= (-(R2+2.5)/(2*L)+unp.sqrt((R2+2.5)**2/(4*L**2)+1/(L*C)))/(2*np.pi)

# Berechnung relativer Fehler nu 1

RF_nu1= (nu1 - nu1_lit.n)/nu1_lit.n

# Berechnung nu 2 starke Dämfung

nu2_lit= ((R2+2.5)/(2*L)+unp.sqrt((R2+2.5)**2/(4*L**2)+1/(L*C)))/(2*np.pi)

# Berechnung relativer Fehler nu 2

RF_nu2= (nu2 - nu2_lit.n)/nu2_lit.n

# Berechnung nu_resonanz

nu_res_lit= (unp.sqrt(1/(L*C)))/(2*np.pi)

#Berechnung relativer Fehler nu_resonanz

RF_nu_res= (nu_res-nu_res_lit.n)/nu_res_lit.n


#Plot für a) #################################################################################

slope, intercept, r_value, p_value, std_err = stats.linregress(ta, np.log(Aa))

#Berechnung T_ex
slopefehl=ufloat(slope, std_err)

T_ex= -1/slopefehl

#Berechnung T_ex_lit

T_ex_lit= 2 * L/(R1+2.5)

#Berechnung Relativer Fehler T_ex

RF_T_ex=(T_ex-T_ex_lit.n)/T_ex_lit.n

#Berechnung R_eff

R_eff=2*L*(-slope)

#Berechnung relativer Fehler R_eff

RF_R_eff=(R_eff.n-R1.n)/R1.n

# Erstellung Plot a)

plt.plot(ta*1e+04,np.log(Aa), 'rx', label="Messdaten")
x_plot = np.linspace(0, 2e-04, 10)
plt.plot(x_plot*1e04,intercept+slope*x_plot, 'k-', label="Lineare Regression")
plt.legend(loc="best")
plt.xlabel('t / $\mu$s')
plt.ylabel(r'$\log (A)  $')
plt.grid()
plt.tight_layout
plt.savefig('build/plota.pdf')
plt.close()

#Plot für c) #######################################################################################################

popU, pcovU = curve_fit(
    U,
    f[0:18],
    A0[0:18],
    sigma=None,
    absolute_sigma=True,
    p0=[1000, 3.5*1e-03, 5*1e-09]
    )

# Erstellung des Plots c) in ln-Darstellung

plt.plot(f, A0, 'rx', label="Messdaten")
x_plot = np.linspace(10*1e03,60*1e03,10000)
plt.plot(x_plot, U(x_plot,*popU),linestyle='-',label='Nichtlineare Regression')
plt.plot(x_plot, U(x_plot,R2.n,L.n, C ),linestyle='--',label='Theoriekurve')
plt.xscale('log')
plt.xlabel(r'$\nu \:/\: $Hz')
plt.ylabel(r'$U_C/U_0  $')
plt.grid()
plt.tight_layout
plt.savefig('build/plotcln.pdf')
plt.close()

# Erstellung des Plots c) in linearer Darstellung

plt.plot(f, A0, 'rx', label="Messdaten")
x_plot = np.linspace(26*1e03,45*1e03,10000)
plt.plot(x_plot, U(x_plot,*popU),linestyle='-',label='Nichtlineare Regression')
plt.plot(x_plot, U(x_plot,R2.n,L.n, C ),linestyle='--',label='Theoriekurve')
plt.xticks( [27500 ,30000,32500, 35000, 37500, 40000, 42500, 45000 ],
            [27.5, 30.0, 32.5, 35.0, 37.5, 40.0, 42.5, 45.0]
)
plt.legend(loc="best")
plt.axhline(Umax/(U0*np.sqrt(2)), linestyle='--') 
plt.axvline(27.8*1e03,linestyle='--', label=r'$\nu_r$') # nu-
plt.axvline(42.7*1e03,linestyle='--') # nu+
plt.xlabel(r'$\nu \:/\: $kHz')
plt.ylabel(r'$U_C/U_0  $')
plt.grid()
plt.tight_layout
plt.savefig('build/plotclin.pdf')
plt.close()

# Plot für d) ###########################################################################################

popphi, pcovphi = curve_fit(
    phi1,
    f[0:18],
    phi[0:18],
    sigma=None,
    absolute_sigma=True,
    p0=[300, 3.5*1e-03, 5*1e-09]
    )

# Erstellung des Plots d) in ln-Darstellung

plt.plot(f*1e-04, phi, 'rx', label="Messdaten")
x_plot = np.linspace(2*1e04,4.5*1e04,10000)
plt.plot(x_plot*1e-04, phi1(x_plot,*popphi),linestyle='-',label='Nichtlineare Regression')
plt.plot(x_plot*1e-04, phi1(x_plot,R2.n, L.n, C),linestyle='-',label='Theoriekurve')
plt.xscale('log')
plt.yticks( [-np.pi/2,-np.pi/4,0,np.pi/4 ,np.pi/2,3 * np.pi/4 ,np.pi],
            [r'$-\pi/2$',r'$-\pi/4$',r'$0$' ,r'$+\pi/4$' ,r'$+\pi/2$' ,r'$+3\pi/4$'  ,r'$+\pi$']
    )
plt.xticks( [3, 4],
            [3, 4]
)
plt.legend(loc="best")
plt.xlabel(r'$\nu \:/\: 10^4$Hz')
plt.ylabel(r'$\phi \:/\: $rad')
plt.tight_layout
plt.grid()
plt.savefig('build/plotdln.pdf')
plt.close()

#Erstellung des Plots d) in linearer Darstellung

plt.plot(f*1e-04, phi, 'rx', label="Messdaten")
x_plot = np.linspace(2*1e04,4.5*1e04,10000)
plt.plot(x_plot*1e-04, phi1(x_plot,*popphi),linestyle='-',label='Nichtlineare Regression')
plt.plot(x_plot*1e-04, phi1(x_plot,R2.n, L.n, C),linestyle='-',label='Theoriekurve')
plt.yticks( [-np.pi/2,-np.pi/4,0,np.pi/4 ,np.pi/2,3 * np.pi/4 ,np.pi],
            [r'$-\pi/2$',r'$-\pi/4$',r'$0$' ,r'$+\pi/4$' ,r'$+\pi/2$' ,r'$+3\pi/4$'  ,r'$+\pi$']
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

# Tex file of R_ex

with open('build/R_eff.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{(R_eff.n):.2f}\pm{R_eff.s:.2f}') 
    RC.write('}{\ohm}$')

# Tex file of RF_R_eff

with open('build/RF_R_eff.tex', 'w') as RC:
    RC.write('$\\num{')
    RC.write(f'{RF_R_eff:.2f}')
    RC.write('}$')

# Tex file of Rap

with open('build/Rap.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{(R*1e-03):.2f}') 
    RC.write('}{\kilo\ohm}$')

# Tex file of nuGraph

with open('build/nuGraph.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{(nuGraph*1e-03):.2f}') 
    RC.write('}{\kilo\hertz}$')

# Tex file of nu1

with open('build/nu1.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{(nu1*1e-03):.2f}') 
    RC.write('}{\kilo\hertz}$')

# Tex file nu2

with open('build/nu2.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{(nu2*1e-03):.2f}') 
    RC.write('}{\kilo\hertz}$')

# Tex file of nu_res

with open('build/nu_res.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{(nu_res*1e-03):.2f}') 
    RC.write('}{\kilo\hertz}$')

# Tex file of L Literaturwert

with open('build/L.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{(L.n*1e03):.2f}\pm{L.s*1e03:.2f}') 
    RC.write('}{\milli\henry}$')

# Tex file of R1 Literaturwert

with open('build/R1.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{(R1.n):.2f}\pm{R1.s:.2f}') 
    RC.write('}{\ohm}$')

# Tex file of R2 Literaturwert

with open('build/R2.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{(R2.n):.2f}\pm{R2.s:.2f}') 
    RC.write('}{\ohm}$')

# Tex file C Literaturwert

with open('build/C.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{C*1e09:.2f}') 
    RC.write('}{\\nano\\farad}$')

# Tex file of R_ap_lit

with open('build/R_ap_lit.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{(R_ap_lit.n*1e-03):.2f}\pm{R_ap_lit.s*1e-03:.2f}') 
    RC.write('}{\kilo\ohm}$')

# Tex file of RF_R

with open('build/RF_R.tex', 'w') as RC:
    RC.write('$\\num{')
    RC.write(f'{RF_R:.4f}')
    RC.write('}$')

# Tex file of nuRech

with open('build/nuRech.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{nuRech.n*1e-03:.2f}\pm{nuRech.s*1e-03:.2f}')
    RC.write('}{\kilo\hertz}$')

# Tex file of q

with open('build/q.tex', 'w') as RC:
    RC.write('$\\num{')
    RC.write(f'{q:.2f}')
    RC.write('}$')

# Tex file of q_lit

with open('build/q_lit.tex', 'w') as RC:
   RC.write('$\\num{')
   RC.write(f'{q_lit.n:.2f}\pm{q_lit.s:.2f}')
   RC.write('}$')

# Tex file of RF_q

with open('build/RF_q.tex', 'w') as RC:
    RC.write('$\\num{')
    RC.write(f'{RF_q:.4f}')
    RC.write('}$')

# Tex file of RF_nu    

with open('build/RF_nu.tex', 'w') as RC:
    RC.write('$\\num{')
    RC.write(f'{RF_nu:.4f}')
    RC.write('}$')

# Tex file of nu1_lit

with open('build/nu1_lit.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{nu1_lit.n*1e-03:.2f}\pm{nu1_lit.s*1e-03:.2f}')
    RC.write('}{\kilo\hertz}$')

# Tex file of RF_nu1

with open('build/RF_nu1.tex', 'w') as RC:
    RC.write('$\\num{')
    RC.write(f'{RF_nu1:.4f}')
    RC.write('}$')

# Tex file of nu2_lit

with open('build/nu2_lit.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{nu2_lit.n*1e-03:.2f}\pm{nu2_lit.s*1e-03:.2f}')
    RC.write('}{\kilo\hertz}$')

# Tex file of RF_nu2

with open('build/RF_nu2.tex', 'w') as RC:
    RC.write('$\\num{')
    RC.write(f'{RF_nu2:.4f}')
    RC.write('}$')

# Tex file of nu_res_lit

with open('build/nu_res_lit.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{nu_res_lit.n*1e-03:.2f}\pm{nu_res_lit.s*1e-03:.2f}')
    RC.write('}{\kilo\hertz}$')

# Tex file of RF_nu_res

with open('build/RF_nu_res.tex', 'w') as RC:
    RC.write('$\\num{')
    RC.write(f'{RF_nu_res:.4f}')
    RC.write('}$')

# Tex file of T_ex

with open('build/T_ex.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{T_ex.n*1e06:.2f}\pm{T_ex.s*1e06:.2f}')
    RC.write('}{\micro\second}$')

with open('build/m.tex', 'w') as RC:
    RC.write('$\\num{')
    RC.write(f'{slope:.2f}\pm{std_err:.2f}')
    RC.write('}$')

with open('build/b.tex', 'w') as RC:
    RC.write('$\\num{')
    RC.write(f'{intercept:.2f}')
    RC.write('}$')

with open('build/T_ex_lit.tex', 'w') as RC:
    RC.write('$\SI{')
    RC.write(f'{T_ex_lit.n*1e06:.2f}\pm{T_ex_lit.s*1e06:.2f}')
    RC.write('}{\micro\second}$')

with open('build/RF_T_ex.tex', 'w') as RC:
    RC.write('$\\num{')
    RC.write(f'{RF_T_ex.n:.4f}')
    RC.write('}$')

# Erstellung Tabelle a) ###################################################################################

Aa1, Aa2 = np.array_split(Aa, 2)
ta1, ta2 = np.array_split(ta*1e+06, 2)

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    {$U_C \:/\: \si{\volt}$} & {$t \:/\: \si{\micro\second}$} &
    {$U_C \:/\: \si{\volt}$} & {$t \:/\: \si{\micro\second}$} \\

    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f} & {2:1.2f} & {3:1.2f}  \\'

# Tabelle für 4a) wird im Tex Format geschrieben ############################################################################################################################################################################

with open('build/table_a.tex', 'w') as g:
    g.write(table_header)
    for row in zip(Aa1, ta1, Aa2, ta2):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

# Erstellung Tabelle c), d) ###################################################################################

f1, f2, f3= np.array_split(f*1e-03, 3)
A1, A2, A3 = np.array_split(A, 3)
t1, t2, t3 = np.array_split(t*1e+06, 3)


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
row_template = r'     {0:1.2f} & {1:1.2f} & {2:1.2f} & {3:1.2f} & {4:1.2f} & {5:1.2f} & {6:1.2f} & {7:1.2f} & {8:1.2f}  \\'

# Tabelle für 4c-d) wird im Tex Format geschrieben ############################################################################################################################################################################

with open('build/table_c.tex', 'w') as h:
    h.write(table_header)
    for row in zip(f1, A1, t1, f2, A2, t2, f3, A3, t3):
        h.write(row_template.format(*row).replace('nan', ''))
        h.write('\n')
    h.write(table_footer)



# Kontrollprints
print(T_ex_lit)
print(1/slope)
print(popU)
print(popphi)

