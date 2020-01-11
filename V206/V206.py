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
        x = '{0:.{1:}f}'.format(num, figures)

    return r'\SI{{{}{}}}{{{}}}'.format(x, exp, unit)

def N(kappa, pb, pa, rho, dmdt):
    return 1/(kappa-1)*(pb*(pa/pb)**(1/kappa)-pa)*1/rho * dmdt


# Definition der Messdaten ##############################################################################################

t, pa, T2, pb, T1, P = np.genfromtxt('data.txt', unpack=True)

t_s = t*60
pa_Pa = pa*1e05+1e05
pb_Pa = pb*1e05+1e05
T2_K = const.convert_temperature(T2, 'c', 'K')
T1_K = const.convert_temperature(T1, 'c', 'K')

P_mean= ufloat(np.mean(P), stats.sem(P))

kappa = 1.14
rho = 5.51 


# Erstellung der Plots ######################################################################################################

# Plot für T1(t) #############################################################################################
par, covm = np.polyfit(t_s, T1_K, deg=2, cov=True)
err = np.sqrt(np.diag(covm))


plt.plot(t_s, T1_K,'kx', label='Messwerte')
x_plot = np.linspace(0, 1800, 5000)
plt.plot(x_plot, par[0]*x_plot**2 + par[1]*x_plot + par[2], 'r-', label="Lineare Regression")
plt.legend(loc="best")
plt.xlabel(r'Zeit $t/s$')
plt.ylabel(r'Temperatur $T_1/K$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotT1.pdf')
plt.close()

parT1=unp.uarray(par, err)

dT1dt = np.array([ 2*parT1[0]*29+parT1[1], 2*parT1[0]*25+parT1[1], 2*parT1[0]*15+parT1[1],  2*parT1[0]*4+parT1[1]])

# Berechnung der aufgenommenen Wärmemenge von T2

dQ1dt=dT1dt * 750*2


# Plot für T2(t) #############################################################################################
par, covm = np.polyfit(t_s, T2_K, deg=2, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(t_s, T2_K,'kx', label='Messwerte')
x_plot = np.linspace(0, 1800, 5000)
plt.plot(x_plot, par[0]*x_plot**2+par[1]*x_plot+par[2], 'r-', label="Lineare Regression")
plt.legend(loc="best")
plt.xlabel(r'Zeit $t/s$')
plt.ylabel(r'Temperatur $T_2/K$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotT2.pdf')
plt.close()

parT2=unp.uarray(par, err)

# Berechnung dT2/dt anhand von 4 Werten

dT2dt = np.array([ 2*parT2[0]*29+parT2[1], 2*parT2[0]*25+parT2[1],2*parT2[0]*15+parT2[1],  2*parT2[0]*4+parT2[1]])

# Berechnung der abgegebenen Wärmemenge von T2

dQ2dt=dT2dt * 750*2


# Plot zur Berechnung von L #####################################################################################

par, covm = np.polyfit(1/(T2_K), np.log(pa_Pa), deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(1/T2_K, np.log(pa_Pa), 'kx', label='Messwerte')
x_plot = np.linspace(0.0033, 0.0037, 1000)
plt.plot(x_plot ,par[0]*x_plot+par[1], 'r-', label="Lineare Regression")
plt.xticks([3.3 *1e-03, 3.4e-03, 3.5e-03, 3.6e-03, 3.7e-03],
           [ 3.3, 3.4, 3.5, 3.6, 3.7])
plt.legend(loc="best")
plt.xlabel(r'$1/T\cdot 10^{-3}\:/\:(1/K) $')
plt.ylabel(r' logarithmischer Dampfdruck $\ln(p/(1\; Pa))$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotL.pdf')
plt.close()

parL=unp.uarray(par, err)
L= -parL[0]* const.R

# weiter Berechnungen ############################################################################################

# Berechnung der reellen Güteziffer

nu_real = dQ1dt/P_mean


# Berecnung der idealen Güteziffer

nu_id = np.array([T1_K[29]/(T1_K[29] -T2_K[29] ), T1_K[25]/(T1_K[25] -T2_K[25] ), T1_K[15]/(T1_K[15] -T2_K[15] ), T1_K[4]/(T1_K[4] -T2_K[4] )])


# Berechnung des Massendurchsatzes

dmdt = 1/L * dQ2dt
print(dmdt)

# Berechnung der mechanischen Leistung des Kompressors

N = np.array([N(kappa, pb_Pa[29], pa_Pa[29], rho, dmdt[0]) ,N(kappa, pb_Pa[25], pa_Pa[25], rho, dmdt[1]) ,N(kappa, pb_Pa[15], pa_Pa[15], rho, dmdt[2]) ,N(kappa, pb_Pa[4], pa_Pa[4], rho, dmdt[3]) ])


# Tex-Dateien ###################################################################################################################

# tex file of parT1_a

with open('build/parT1_a.tex', 'w') as f: 
  f.write(make_SI(parT1[0],r'' ,figures=1))

# tex file of parT1_b

with open('build/parT1_b.tex', 'w') as f: 
  f.write(make_SI(parT1[1],r'' ,figures=1))

# tex file of parT1_c

with open('build/parT1_c.tex', 'w') as f: 
  f.write(make_SI(parT1[2],r'' ,figures=1))

# tex file of parT2_a

with open('build/parT2_a.tex', 'w') as f: 
  f.write(make_SI(parT1[0],r'' ,figures=1))

# tex file of parT1_b

with open('build/parT2_b.tex', 'w') as f: 
  f.write(make_SI(parT1[1],r'' ,figures=1))

# tex file of parT2_c

with open('build/parT2_c.tex', 'w') as f: 
  f.write(make_SI(parT1[2],r'' ,figures=1))

# tex file of parL_m

with open('build/parL_m.tex', 'w') as f: 
  f.write(make_SI(parL[0],r'' ,figures=1))

# tex file of L 

with open('build/parL_b.tex', 'w') as f: 
  f.write(make_SI(parL[1],r'' ,figures=1))


# tex file of L 

with open('build/L.tex', 'w') as f: 
  f.write(make_SI(L,r'\kilo\joule\mol\tothe{-1}' ,figures=1))


# Tabellen ########################################################################################################################################################

# Tabelle für calc.tex -------------------------------------------------------------------------------------------------------------------------------------------------

t_rech = np.array([29, 25, 15, 4])

table_header = r'''
  \begin{tabular}{c c c c c c}
    \toprule
     \multicolumn{1}{c}{Zeit t$\:/\: \si{\second}$} &\multicolumn{1}{c}{$\frac{dQ_1}{dt}$} & \multicolumn{1}{c}{$\nu_{real}\cdot 10^ {-3} $} & \multicolumn{1}{c}{$\nu_{ideal}$} & \multicolumn{1}{c}{$\frac{dQ_2}{dt}$ } & \multicolumn{1}{c}{$\frac{dm}{dt}\cdot 10^{-4}$ }\\
      
    \cmidrule(lr){1-6}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.0f} & {1:1.1f} & {2:1.1f} & {3:1.1f} & {4:1.1f} & {5:1.2f} \\'


with open('build/table_calc.tex', 'w') as g:
    g.write(table_header)
    for row in zip(t_rech, dQ1dt, nu_real*1e03, nu_id, dQ2dt, dmdt*1e04):
        g.write(row_template.format(*row).replace('nan', '').replace('+/-','\pm'))
        g.write('\n')
    g.write(table_footer)

# Tabelle der Messwerte

table_header = r'''
  \begin{tabular}{c c c c c c}
    \toprule
     \multicolumn{1}{c}{Zeit} &\multicolumn{1}{c}{Druck} & \multicolumn{1}{c}{Temperatur} & \multicolumn{1}{c}{Temperatur} & \multicolumn{1}{c}{Druck} & \multicolumn{1}{c}{Leistung }\\
     \multicolumn{1}{c}{$t\:/\: \si{\second}$} &\multicolumn{1}{c}{$p_a \:/\: \si{\bar} $} & \multicolumn{1}{c}{$ T_2 \:/\: \si{\celsius} $} & \multicolumn{1}{c}{$T_1 \:/\; \si{\celsius} $} & \multicolumn{1}{c}{$p_b\:/\: \si{\bar} $ } & \multicolumn{1}{c}{$P\:/\: \si{\watt} $ }\\

    \cmidrule(lr){1-6}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.0f} & {1:1.1f} & {2:1.1f} & {3:1.1f} & {4:1.1f} & {5:1.2f} \\'


with open('build/table.tex', 'w') as g:
    g.write(table_header)
    for row in zip(t, pa, T2, pb, T1, P):
        g.write(row_template.format(*row).replace('nan', ''))
        g.write('\n')
    g.write(table_footer)

