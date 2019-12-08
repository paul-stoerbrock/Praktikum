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

def I(T, D, I_D, m, a):
  return (T/(2*np.pi))**2 * D - I_D - m * a**2

def make_SI(num, unit, exp='', figures=None):
    ''' Format an uncertainties ufloat as a \SI quantity '''
    if np.any(stds([num])):
        if figures is None:
            figures = ''
        x = '{0:.{1:}uf}'.format(num, figures).replace('/', '')
    else:
        x = '{0:.{1:}f}'.format(num, figures)

    return r'\SI{{{}{}}}{{{}}}'.format(x, exp, unit)

# Messwerte

F, phi_Grad, a, T_I_Stab1, T_I_Stab2, T_Kugel, T_Zylinder, T_Puppe_fast, T_Puppe_slow = np.genfromtxt('data.txt', unpack=True) 
# F=Federkraft, phi_Grad=Winkel in Gradmaß, a=Abstand Drehachse(rot,blau), 
#T_I_Stab1 und 2 = dreifache Periodendauer(zur Bestimmung Trägheitsmoment Stab),
#T_Kugel = Periodendauer Kugel, T_Zylinder =  Periodendauer des Zylinders, 
#T_Puppe_fast =  Periodendauer Puppe angelegte Arme,
#T_Puppe_slow = Periodendauer Puppe ausgestreckte Arme

# Berechnungen ##########################################################################################################

phi_Bogen = phi_Grad/180 * np.pi # Umrechnung von Gradmaß in Bogenmaß

T_I_Stab = (np.round(T_I_Stab1,1) + np.round(T_I_Stab2))/6

D = (F*30)/phi_Bogen # Bestimmung der Winkelrichtgröße

D_mittelw = np.mean(D) # Bestimmung der Federkonstante



# Lineare Regression zur Bestimmung des Trägheitsmoments der Drehachse

slope, intercept, r_value, p_value, std_err = stats.linregress((a[0:10])**2 , T_I_Stab[0:10]**2)



# Berechnung der Trägheitsmomente #####################################################################

# experimenteller Wert des Trägheitsmoment der Drillachse

I_D = intercept * (4*np.pi**2)/D_mittelw *1e-05 

# experimenteller Wert des Trägheitsmoment des Zylinders

I_Zylinder = I(np.round(T_Zylinder[0:5],1), D_mittelw, I_D, 0, 0)
I_Zylinder_mean = np.mean(I_Zylinder)

# theoretischer Wert des Trägheitsmoments des Zylinders

I_Zylinder_Theorie = 1525.5*((4**2)/4+(13.95**2)/12)*1e-05

# experimenteller Wert des Trägheitsmoment der Kugel

I_Kugel = I(np.round(T_Kugel[0:5],1), D_mittelw, I_D, 0, 0)
I_Kugel_mean = np.mean(I_Kugel)

# theoretischer Wert des Trägheitsmoments der Kugel

I_Kugel_Theorie = 2/5 * 1168.5 * (14.6/2)**2 *1e-05

# theoretischer Wert des Trägheitsmoment der Puppe Arme angelegt

m_Torso = np.pi * (3.89/2)**2 * 9.625 * 0.78
m_Arm  = np.pi * (1.65/2)**2 * 13.25 * 0.78
m_Kopf = np.pi * (3.265/2)**2 * 4.380 * 0.78
m_Bein = np.pi * (2.18/2)**2 * 14.475 * 0.78

I_Puppe_an_theo = (m_Torso * (3.89/2)**2/2 + m_Kopf * (3.265/2)**2/2 + 2 * m_Bein * (2.18/2)**2/2 + 2 * ( m_Arm * (1.65/2)**2/2 + m_Arm * (1.65/2 + 3.89/2)**2))*1e-05

# experimenteller Wert des Trägheitsmoment der Puppe Arme angelegt

I_Puppe_an_exp = (2 * I(T_Puppe_fast[0:5], D_mittelw, I_D, m_Arm, 1.65/2 + 3.89/2) + I(T_Puppe_fast[0:5], D_mittelw, I_D, m_Kopf, 0) + I(T_Puppe_fast[0:5], D_mittelw, I_D, m_Torso, 0) + 2* I(T_Puppe_fast[0:5], D_mittelw, I_D,2 * m_Bein, 2.18/2))*1e-05
I_Puppe_an_exp_mean = np.mean(I_Puppe_an_exp)

# theoretischer Wert des Trägheitsmoment der Puppe Arme ausgestreckt

I_Puppe_aus_theo = (m_Torso * (3.89/2)**2/2 + m_Kopf * (3.265/2)**2/2 + 2 * m_Bein * (2.18/2)**2/2 + 2 * ( m_Arm * ((1.65/2)**2/4+(13.25)**2/12) + m_Arm * (13.25/2 + 3.89/2)**2))*1e-05

# experimanteller Wert des Trägheitsmoment der Puppe Arme ausgestreckt

I_Puppe_aus_exp =(2 * I(T_Puppe_slow[0:5], D_mittelw, I_D, m_Arm, 13.25/2 + 3.89/2) + I(T_Puppe_slow[0:5], D_mittelw, I_D, m_Kopf, 0) + I(T_Puppe_slow[0:5], D_mittelw, I_D, m_Torso, 0) + 2* (I(T_Puppe_slow[0:5], D_mittelw, I_D,2 * m_Bein, 2.18/2)))*1e-05
I_Puppe_aus_exp_mean = np.mean(I_Puppe_aus_exp)


# Berechnung relativer Fehler #####################################################################################################################################

# Relativer Fehler von I_Zylinder

RF_I_Zylinder = (I_Zylinder_mean - I_Zylinder_Theorie)/I_Zylinder_Theorie

# Relativer Fehler von I_Kugel

RF_I_Kugel = (I_Kugel_mean-I_Kugel_Theorie)/I_Kugel_Theorie

# Relativer Fehler von I_Puppe_an

RF_I_Puppe_an = (I_Puppe_an_exp_mean-I_Puppe_an_theo)/I_Puppe_an_theo

# Relativer Fehler von I_Puppe_aus

RF_I_Puppe_aus = (I_Puppe_aus_exp_mean-I_Puppe_aus_theo)/I_Puppe_aus_theo


# Erstellung der Plots ######################################################################################################################

# Plot zur Bestimmung des Trägheitsmoment I_Stab

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

# SI-tex-files ####################################################################################################################

# tex file of D_mittelw

with open('build/D_mean.tex', 'w') as f:
  f.write(make_SI(D_mittelw , r'\newton\meter', figures=2))

# tex file of slope

with open ('build/slope.tex', 'w') as f:
  f.write(make_SI(slope, '', figures=2))

# tex file of intercept

with open('build/intercept.tex', 'w') as f:
  f.write(make_SI(intercept, '', figures=2))

# tex file of I_D

with open ('build/I_D.tex', 'w') as f:
  f.write(make_SI(I_D*1e03, r'\kilo\gram\meter\square',  exp='1e-03', figures=2))

# tex file of I_Zylinder_mean

with open ('build/I_Zylinder_mean.tex', 'w') as f:
  f.write(make_SI(I_Zylinder_mean, r'\kilo\gram\meter\square', figures=2))

# tex file of I_zylinder_Theorie

with open ('build/I_Zylinder_Theorie.tex', 'w') as f:
  f.write(make_SI(I_Zylinder_Theorie, r'\kilo\gram\meter\square', figures=2))

# tex file of I_Kugel_mean

with open ('build/I_Kugel_mean.tex', 'w') as f:
  f.write(make_SI(I_Kugel_mean, r'\kilo\gram\meter\square', figures=2))

# tex file of I_Kugel_Theorie

with open ('build/I_Kugel_Theorie.tex', 'w') as f:
  f.write(make_SI(I_Kugel_Theorie, r'\kilo\gram\meter\square', figures=2))

#texfile of m_Torso

with open ('build/m_Torso.tex', 'w') as f:
  f.write(make_SI(m_Torso, r'\gram', figures=2))

# tex file of m_Arm

with open ('build/m_Arm.tex', 'w') as f:
  f.write(make_SI(m_Arm, r'\gram', figures=2))

# tex file of m_Kopf

with open ('build/m_Kopf.tex', 'w') as f:
  f.write(make_SI(m_Kopf, r'\gram', figures=2))

# tex file of m_Bein

with open ('build/m_Bein.tex', 'w') as f:
  f.write(make_SI(m_Bein, r'\gram', figures=2))

# tex file of I_Puppe_an_theo

with open ('build/I_Puppe_an_theo.tex', 'w') as f:
  f.write(make_SI(I_Puppe_an_theo*1e03, r'\kilo\gram\meter\square', exp='1e-03', figures=1))

# tex file of I_Puppe_an_exp_mean

with open ('build/I_Puppe_an_exp_mean.tex', 'w') as f:
  f.write(make_SI(I_Puppe_an_exp_mean*1e03, r'\kilo\gram\meter\square', exp='1e-03',figures=1))

# tex file of I_Puppe_aus_theo

with open ('build/I_Puppe_aus_theo.tex', 'w') as f:
  f.write(make_SI(I_Puppe_aus_theo*1e03, r'\kilo\gram\meter\square', exp='1e-03',figures=1))

# tex file of I_Puppe_aus_exp_mean

with open ('build/I_Puppe_aus_exp_mean.tex', 'w') as f:
  f.write(make_SI(I_Puppe_aus_exp_mean*1e03, r'\kilo\gram\meter\square', exp='1e-03',figures=1))

# tex file of RF_I_Zylinder

with open ('build/RF_I_Zylinder.tex', 'w') as f:
  f.write(make_SI(RF_I_Zylinder, r'', figures=2))

# tex file of RF_I_Kugel

with open ('build/RF_I_Kugel.tex', 'w') as f:
  f.write(make_SI(RF_I_Kugel, r'', figures=2))

# tex file of RF_I_Puppe_an

with open ('build/RF_I_Puppe_an.tex', 'w') as f:
  f.write(make_SI(RF_I_Puppe_an, r'', figures=2))

# tex file of RF_I_Puppe_aus

with open ('build/RF_I_Puppe_aus.tex', 'w') as f:
  f.write(make_SI(RF_I_Puppe_aus, r'', figures=2))


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
row_template = r'     {0:1.0f} & {1:1.2f} & {2:1.0f} & {3:1.2f}  \\'

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
  \begin{tabular}{c c c c c c}
    \toprule
    {$a \:/\: \si{\centi\meter}$} & {$T_{Wert1} \:/\: \si{\second}$} & {$T_{Wert2} \:/\: \si{\second}$} &
    {$a \:/\: \si{\centi\meter}$} & {$T_{Wert1} \:/\: \si{\second}$} & {$T_{Wert2} \:/\: \si{\second}$} \\
   

    \cmidrule(lr{0,5em}){1-3} \cmidrule(lr{0,5em}){4-6}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f} & {2:1.2f} & {3:1.2f} & {4:1.2f} & {5:1.2f}  \\'

# Tabelle für Trägheitsmoment Stab wird im Tex Format geschrieben ############################################################################################################################################################################
a1,a2 = np.array_split(a[0:10],2)
T_I_Stab11, T_I_Stab12 = np.array_split(T_I_Stab1[0:10], 2)
T_I_Stab21, T_I_Stab22 = np.array_split(T_I_Stab2[0:10], 2)

with open('build/table_IStab.tex', 'w') as g:
    g.write(table_header)
    for row in zip(a1, T_I_Stab11,T_I_Stab21, a2, T_I_Stab12, T_I_Stab22):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

#Tabelle für Trägheitsmoment Kugel, Zylinder und Puppe

table_header = r'''
  \begin{tabular}{c c c c}
    \toprule
    {$T_{\text{Kugel}} \:/\: \si{\second}$} & {$T_{\text{Zylinder}} \:/\: \si{\second}$} &
    {$T_{\text{Puppe,an}} \:/\: \si{\second}$} & {$T_{\text{Puppe,aus}} \:/\: \si{\second}$} \\

    \cmidrule(lr{0,5em}){1-1} \cmidrule(lr{0,5em}){2-2} \cmidrule(lr{0,5em}){3-3} \cmidrule(lr{0,5em}){4-4}
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

print(I_Kugel_mean)
print(I_Kugel_Theorie)
print(I_Zylinder_mean)
print(I_Zylinder_Theorie)



