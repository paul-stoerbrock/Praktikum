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

# Messwerte

F, phi_Grad, a, T_I_Stab1, T_I_Stab2, T_Kugel, T_Zylinder, T_Puppe_fast, T_Puppe_slow = np.genfromtxt('data.txt', unpack=True) 
# F=Federkraft, phi_Grad=Winkel in Gradmaß, a=Abstand Drehachse(rot,blau), 
#T_I_Stab1 und 2 = dreifache Periodendauer(zu Bestimmun Trägheitsmoment Stab),
#T_Kugel = Periodendauer Kugel, T_Zylinder =  Periodendauer des Zylinders, 
#T_Puppe_fast =  Periodendauer Puppe angelegte Arme,
#T_Puppe_slow = Periodendauer Puppe ausgestreckte Arme

phi_Bogen = phi_Grad/180 * np.pi # Umrechnung von Gradmaß in Bogenmaß

T_I_Stab = (T_I_Stab1 + T_I_Stab2)/6

# Erstellung der Plots ######################################################################################################################

# Plot zur Bestimmung des Trägheitsmoment I_Stab

slope, intercept, r_value, p_value, std_err = stats.linregress(a[0:10]**2 , T_I_Stab[0:10]**2)

plt.plot(a[0:10]**2 , T_I_Stab[0:10]**2, 'bx', label="Messdaten")
x_plot = np.linspace(0, 700, 1000)
plt.plot(x_plot,intercept+slope*x_plot, 'k-', label="Lineare Regression")
plt.legend(loc="best")
plt.xlabel(r'$a^2$/${cm}^2$')
plt.ylabel(r'$T^2$/$s^2$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot.pdf')
plt.close()

# Tabellen #############################################################################################################



#Tabelle für Winkelrichtgröße D

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    {$\varphi \:/\: \text{grad}$} & {$F\:/\: \si{\newton}$} &
    {$\varphi \:/\: \text{grad}$} & {$F\:/\: \si{\newton}$} \\

    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f} & {2:1.2f} & {3:1.2f}  \\'

# Tabelle für Winkelrichtgröße D wird im Tex Format geschrieben ############################################################################################################################################################################

F1, F2 = np.array_split(F, 2)
phi1, phi2 = np.array_split(phi_Grad, 2)

with open('build/table_D.tex', 'w') as g:
    g.write(table_header)
    for row in zip(phi1, F1, phi2, F2):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)


#Tabelle für Trägheitsmoment Stab

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    {$a \:/\: \si{\centi\meter}$} & {$T_{Stab} \:/\: \si{\second}$} &
    {$a \:/\: \si{\centi\meter}$} & {$T_{Stab} \:/\: \si{\second}$} \\
   

    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f} & {2:1.2f} & {3:1.2f}  \\'

# Tabelle für Trägheitsmoment Stab wird im Tex Format geschrieben ############################################################################################################################################################################
a1,a2 = np.array_split(a[0:10],2)
T_I_Stab1, T_I_Stab2 = np.array_split(T_I_Stab[0:10], 2)

with open('build/table_IStab.tex', 'w') as g:
    g.write(table_header)
    for row in zip(a1, T_I_Stab1, a2, T_I_Stab2):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

#Tabelle für Trägheitsmoment Kugel, Zylinder und Puppe

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    {$T_{\text{Kugel}} \:/\: \si{\second}$} & {$T_{\text{Zylinder}} \:/\: \si{\second}$} &
    {$T_{\text{Puppe,an}} \:/\: \si{\second}$} & {$T_{\text{Puppe,aus}} \:/\: \si{\second}$} \\

    \cmidrule(lr{0,5em}){1} \cmidrule(lr{0,5em}){2} \cmidrule(lr{0,5em}){3} \cmidrule(lr{0,5em}){4}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f} & {2:1.2f} & {3:1.2f}  \\'

# Tabelle für Trägheitsmoment Kugel, Zylinder und Puppe wird im Tex Format geschrieben ############################################################################################################################################################################


with open('build/table_I.tex', 'w') as g:
    g.write(table_header)
    for row in zip(T_Kugel[0:5], T_Zylinder[0:5], T_Puppe_fast[0:5], T_Puppe_slow[0:5]):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

# Testprints ##########################################################################################

