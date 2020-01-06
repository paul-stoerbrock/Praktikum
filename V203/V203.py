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

# Definition der Variablen ########################################################################################################################################################

p_k1bar, C_k1bar = np.genfromtxt('datak1b.txt', unpack=True)
p_g1bar, C_g1bar = np.genfromtxt('datag1b.txt', unpack=True)

T = 373 # Temperatur für größer 1 bar

K_k1bar = const.convert_temperature(C_k1bar[:77], 'c', 'K')

K_g1bar = const.convert_temperature(C_g1bar[:77], 'c', 'K') 

pinP_k1bar = p_k1bar[:77] * 1e05
pinP_g1bar = p_g1bar[:77] * 1e05 

# Erstellung der Plots #################################################################################################################################################################

# Plot für kleiner 1 bar

park1b, covmk1b = np.polyfit(1/(K_k1bar), np.log(pinP_k1bar[:77]), deg=1, cov=True)
errk1b = np.sqrt(np.diag(covmk1b))

plt.plot(1/K_k1bar, np.log(pinP_k1bar[:77]), 'kx', label='Messwerte')
x_plot = np.linspace(0.0025, 0.0035, 1000)
plt.plot(x_plot,park1b[1] +park1b[0]*x_plot, 'r-', label="Lineare Regression")
plt.xticks([2.6*1e-03, 2.8*1e-03, 3*1e-03, 3.2*1e-03,3.4 *1e-03],
           [2.6, 2.8, 3, 3.2, 3.4])
plt.legend(loc="best")
plt.xlabel(r'$1/T\cdot 10^{-3}\:/\:(1/K) $')
plt.ylabel(r' logarithmischer Dampfdruck $\ln(p/(1\; Pa))$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotk1b.pdf')
plt.close()

park1b=unp.uarray(park1b, errk1b)
Lk1b= -park1b[0]* const.R

L_a = const.R * T
L_i = (Lk1b-L_a)*0.01036


# Plot für größer 1 bar ############################################################################################

parg1b, covmg1b = np.polyfit(K_g1bar, pinP_g1bar, deg=3, cov=True)
err = np.sqrt(np.diag(covmg1b))

# Testplot, um zu sehen ob Kurve gut ist

plt.plot(K_g1bar, pinP_g1bar, 'kx', label='Messwerte')
x_plot = np.linspace(90, 400, 1000)
plt.plot(x_plot ,parg1b[0]*x_plot**3+parg1b[1]*x_plot**2+parg1b[2]*x_plot*x_plot+parg1b[3] ,'r-' , label="Lineare Regression")
plt.savefig('build/test.pdf')
plt.close()


parg1b = unp.uarray(parg1b ,err)

# Berechnung des Wertes dpdT

dpdT=3*parg1b[0] * K_g1bar**2 + 2*parg1b[1]*K_g1bar+parg1b[2]

# Berechnung des Dampfvolumens

VDplus = (const.R * K_g1bar)/(2*pinP_g1bar) + np.sqrt((const.R**2*K_g1bar**2)/(4*pinP_g1bar**2)-0.9/pinP_g1bar)
VDminus = (const.R * K_g1bar)/(2*pinP_g1bar) - np.sqrt((const.R**2*K_g1bar**2)/(4*pinP_g1bar**2)-0.9/pinP_g1bar)

# Berechnung der Verdampfungswärme

# Für VDplus
Lg1b = dpdT*VDplus*K_g1bar


parg1bL, covmg1b = np.polyfit(K_g1bar, noms(Lg1b), deg=1, cov=True)
err = np.sqrt(np.diag(covmg1b))


plt.plot(K_g1bar, noms(Lg1b), 'kx', label='Messwerte')
x_plot = np.linspace(360, 480, 1000)
plt.plot(x_plot ,parg1bL[0]*x_plot+parg1bL[1] ,'r-' , label="Lineare Regression")
plt.yticks([35*1e03, 40*1e03, 45*1e03, 50*1e03, 55*1e03, 60*1e03],
           [35, 40, 45, 50, 45, 50, 55, 60])
plt.legend(loc="best")
plt.xlabel(r'Temperatur $T\:/\:K$')
plt.ylabel(r'Verdampfungswärme $L\:/\:kJ$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotg1bplus.pdf')
plt.close()

parg1bLplus = unp.uarray(parg1bL ,err)

# Für VDminus
Lg1b = dpdT*VDminus*K_g1bar


parg1bL, covmg1bL = np.polyfit(K_g1bar, noms(Lg1b), deg=1, cov=True)
err = np.sqrt(np.diag(covmg1b))


plt.plot(K_g1bar, noms(Lg1b), 'kx', label='Messwerte')
x_plot = np.linspace(360, 480, 1000)
plt.plot(x_plot ,parg1bL[0]*x_plot+parg1bL[1] ,'r-' , label="Lineare Regression")
#plt.yticks([35*1e03, 40*1e03, 45*1e03, 50*1e03, 55*1e03, 60*1e03],
#           [35, 40, 45, 50, 45, 50, 55, 60])
plt.legend(loc="best")
plt.xlabel(r'Temperatur $T\:/\:K$')
plt.ylabel(r'Verdampfungswärme $L\:/\:kJ$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotg1bminus.pdf')
plt.close()

parg1bLminus = unp.uarray(parg1bL ,err)

# Parameter werden ins tex-Format geschrieben ########################################################################################################################################################

# tex file of m_plotk1b ################################################################################

with open('build/m_plotk1b.tex', 'w') as f: 
  f.write(make_SI(park1b[0]*1e-03,r'' ,exp='e03' ,figures=1))

# tex file of b_plotk1b ###############################################################

with open('build/b_plotk1b.tex', 'w') as f: 
  f.write(make_SI(park1b[1] ,r'', figures=2))

#tex file of Lk1b ##################################################################

with open('build/Lk1b.tex', 'w') as f: 
  f.write(make_SI(Lk1b*1e-03 ,r'\kilo\joule\mol\tothe{-1}', figures=1))

#tex file of L_a ###########################################################################

with open('build/L_a.tex', 'w') as f: 
  f.write(make_SI(L_a*1e-03 ,r'\kilo\joule\mol\tothe{-1}', figures=2))

#tex file of L_i #####################################################################

with open('build/L_i.tex', 'w') as f: 
  f.write(make_SI(L_i*1e-03 ,r'\electronvolt', figures=2))

#tex file of factor a ###########################################################################

with open('build/parg1b_a.tex', 'w') as f: 
  f.write(make_SI(parg1b[0] ,r'', figures=2))

#tex file of factor b ###########################################################################

with open('build/parg1b_b.tex', 'w') as f: 
  f.write(make_SI(parg1b[1] ,r'', figures=2))

#tex file of factor c ###########################################################################

with open('build/parg1b_c.tex', 'w') as f: 
  f.write(make_SI(parg1b[2] ,r'', figures=2))

#tex file of factor d ###########################################################################

with open('build/parg1b_d.tex', 'w') as f: 
  f.write(make_SI(parg1b[3] ,r'', figures=2))

#tex file of dp/dT  ###########################################################################

with open('build/dpdT.tex', 'w') as f: 
  f.write(make_SI(np.mean(dpdT) ,r'', figures=2))

#tex file of VDplus ###########################################################################

with open('build/VDplus.tex', 'w') as f: 
  f.write(make_SI(np.mean(VDplus) ,r'\meter\tothe{3}', figures=2))

#tex file of VDminus ###########################################################################

with open('build/VDminus.tex', 'w') as f: 
  f.write(make_SI(np.mean(VDminus) ,r'\meter\tothe{3}', figures=2))

#tex file of m L_g1bplus ###########################################################################

with open('build/Lg1b_mplus.tex', 'w') as f: 
  f.write(make_SI(parg1bLplus[0] ,r'', figures=2))

#tex file of b L_g1bplus  ###########################################################################

with open('build/Lg1b_bplus.tex', 'w') as f: 
  f.write(make_SI(parg1bLplus[1]*1e-05 ,r'',exp='e05' ,figures=1))

#tex file of m L_g1bminus ###########################################################################

with open('build/Lg1b_mminus.tex', 'w') as f: 
  f.write(make_SI(parg1bLminus[0] ,r'', figures=2))

#tex file of b L_g1bminus  ###########################################################################

with open('build/Lg1b_bminus.tex', 'w') as f: 
  f.write(make_SI(parg1bLminus[1]*1e-05 ,r'',exp='e05' ,figures=1))

# Tabellen ########################################################################################################################################################

# Tabelle für datak1b.tex -------------------------------------------------------------------------------------------------------------------------------------------------

p_k1bar1, p_k1bar2, p_k1bar3 = np.array_split(p_k1bar*1e+03, 3)
C_k1bar1, C_k1bar2, C_k1bar3 = np.array_split(C_k1bar, 3)

table_header = r'''
  \begin{tabular}{c c c c c c}
    \toprule
    \multicolumn{1}{c}{Druck\:/\:} & \multicolumn{1}{c}{Temperatur\:/\:} & \multicolumn{1}{c}{Druck\:/\:} & \multicolumn{1}{c}{Temperatur\:/\:} & \multicolumn{1}{c}{Druck\:/\:} & \multicolumn{1}{c}{Temperatur\:/\:}\\
    \multicolumn{1}{c}{$\si{\milli\bar}$} & \multicolumn{1}{c}{$\si{\celsius}$} & \multicolumn{1}{c}{$\si{\milli\bar}$} & \multicolumn{1}{c}{$\si{\celsius}$} & \multicolumn{1}{c}{$\si{\milli\bar}$} & \multicolumn{1}{c}{$\si{\celsius}$}\\

    \cmidrule(lr){1-2} \cmidrule(lr){3-4} \cmidrule(lr){5-6}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.0f} & {1:1.0f} & {2:1.0f} & {3:1.0f} & {4:1.0f} & {5:1.0f} \\'


with open('build/table_k1bar.tex', 'w') as g:
    g.write(table_header)
    for row in zip(p_k1bar1, C_k1bar1, p_k1bar2, C_k1bar2, p_k1bar3, C_k1bar3):
        g.write(row_template.format(*row).replace('nan', ''))
        g.write('\n')
    g.write(table_footer)

# Tabelle für datag1b.tex -------------------------------------------------------------------------------------------------------------------------------------------------

p_g1bar1, p_g1bar2, p_g1bar3 = np.array_split(p_g1bar, 3)
C_g1bar1, C_g1bar2, C_g1bar3 = np.array_split(C_g1bar, 3)

table_header = r'''
  \begin{tabular}{c c c c c c}
    \toprule
    \multicolumn{1}{c}{Druck\:/\:} & \multicolumn{1}{c}{Temperatur\:/\:} & \multicolumn{1}{c}{Druck\:/\:} & \multicolumn{1}{c}{Temperatur\:/\:} & \multicolumn{1}{c}{Druck\:/\:} & \multicolumn{1}{c}{Temperatur\:/\:}\\
    \multicolumn{1}{c}{$\si{\milli\bar}$} & \multicolumn{1}{c}{$\si{\celsius}$} & \multicolumn{1}{c}{$\si{\milli\bar}$} & \multicolumn{1}{c}{$\si{\celsius}$} & \multicolumn{1}{c}{$\si{\milli\bar}$} & \multicolumn{1}{c}{$\si{\celsius}$}\\
    \cmidrule(lr){1-2} \cmidrule(lr{0,5em}){3-4} \cmidrule(lr{0,5em}){5-6}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.0f} & {2:1.1f} & {3:1.0f} & {4:1.1f} & {5:1.0f} \\'


with open('build/table_g1bar.tex', 'w') as g:
    g.write(table_header)
    for row in zip(p_g1bar1, C_g1bar1, p_g1bar2, C_g1bar2, p_g1bar3, C_g1bar3):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)