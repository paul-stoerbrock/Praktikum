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

def make_SI(num, unit, exp='', figures=None):
    ''' Format an uncertainties ufloat as a \SI quantity '''
    if np.any(stds([num])):
        if figures is None:
            figures = ''
        x = '{0:.{1:}uf}'.format(num, figures).replace('/', '')
    else:
        x = '{0:.{1:}f}'.format(num, figures)

    return r'\SI{{{}{}}}{{{}}}'.format(x, exp, unit)

def x_dopl(x,l):
  return 3*l**2*x-4*x**3

def x_dopr(x,l):
  return -12*l*x**2 + 4*x**3 +9*l**2 * x - l**3

# Messwerte ####################################################################################################################

xCu_einohne, DCu_einohne,xCu_einmit, DCu_einmit = np.genfromtxt('dataCuein.txt', unpack=True) # Messwette für Kupfer einseitig belastet 
xCu_dopohne, DCu_dopohne,xCu_dopmit, DCu_dopmit = np.genfromtxt('dataCudop.txt', unpack=True) # Messwerte für Kupfer doppelseitig belastet
xAl_einohne, DAl_einohne,xAl_einmit, DAl_einmit = np.genfromtxt('dataAlein.txt', unpack=True) # Messwerte für Aluminium einseitig belastet
xAl_dopohne, DAl_dopohne,xAl_dopmit, DAl_dopmit = np.genfromtxt('dataAldop.txt', unpack=True) # Messwerte für Aluminium doppelseitig belastet

xCu_einohne *= 1e-02
xCu_dopohne *= 1e-02
xAl_einohne *= 1e-02
xAl_dopohne *= 1e-02

#Maße der Stäbe
l_Cu = 0.600
l_Al = 0.592
d_Al =np.array([10.00, 10.00, 10.10, 10.00, 10.10, 10.10, 10.10, 10.00, 10.25, 10.00])*1e-03
d_Cu =np.array([10.02, 10.05, 10.06, 10.03, 10.04, 10.09, 10.04, 10.05, 10.08, 10.05])*1e-03
b_Cu =np.array([10.06, 10.02, 10.05, 10.01, 10.03, 10.03, 10.04, 10.02, 10.04, 10.03])*1e-03

# Massen der Körper
m_aufhaeng = 19 *1e-03
m_schraube = 22.1*1e-03
m_Cu_stange = 528.8*1e-03
m_Al_stange = 132.6*1e-03

# Gewichte Cuein
m_Cuein1 = 502.5*1e-03
m_Cuein2 = 503.3*1e-03

# Gewichte Cudop
m_Cudop1 = 1170.5*1e-03
m_Cudop2 = 1159.7*1e-03
m_Cudop3 = 1162.3*1e-03
m_Cudop4 = 1170.8*1e-03

# Gewichte Alein
m_Alein1 = 500.1*1e-03

# Gewichte Aldop
m_Aldop1 = 500.1*1e-03
m_Aldop2 = 499.8*1e-03
m_Aldop3 = 226.7*1e-03

# Berechnungen ###########################################################################################################

# Berechnung der Differenz von DCu_einmit - DCu_einohne

D_Cu_einDiff = (DCu_einohne - DCu_einmit ) *1e-03

# Berechnung der Differenz von DCu_dopmit - DCu_dopohne

D_Cu_dopDiff = (DCu_dopohne - DCu_dopmit ) *1e-03

# Berechnung der Differenz von DAl_einmit - DAl_einohne

D_Al_einDiff = (DAl_einohne - DAl_einmit ) *1e-03

# Berechnung der Differenz von DAl_dopmit - DAl_einohne

D_Al_dopDiff = (DAl_dopohne - DAl_dopmit ) *1e-03

# Berechnung des Flächenträgheitsmoment

# Für Cu

I_Cu = (np.mean(b_Cu)**3*np.mean(d_Cu) )/12

# Für Al

I_Al = (np.mean(d_Al)**4*np.pi)/64


# Erstellung der Plots ###############################################################################################################

# Plot für Kupfer einseitig belastet ###############################################################################################################

# Für Cuein

slopeCuein, interceptCuein, r_valueCuein, p_valueCuein, std_errCuein = stats.linregress(l_Cu * xCu_einohne**2 -xCu_einohne**3/3 , D_Cu_einDiff)

plt.plot(l_Cu * xCu_einohne**2 -xCu_einohne**3/3 ,D_Cu_einDiff , 'bx', label="Messdaten")
x_plot = np.linspace(0, 0.1, 1000)
plt.plot(x_plot,interceptCuein+slopeCuein*x_plot, 'k-', label="Lineare Regression")
plt.yticks( [0 ,1e-03,2e-03, 3e-03, 4e-03, 5e-03],
            [0, 1, 2, 3, 4, 5]
)
plt.legend(loc="best")
plt.xlabel(r'$L-x^2-x^3/3$')
plt.ylabel(r'Durchbiegung D/mm')
plt.grid()
plt.tight_layout
plt.savefig('build/plotCuein.pdf')
plt.close()

#Elastizitätsmodul für Cuein

E_Cuein = ((m_aufhaeng+m_schraube+m_Cuein1+m_Cuein2)*const.g)/(2*I_Cu*slopeCuein)

# Plot für Aluminium einseitig belastet ####################################################################################################################

slopeAlein, interceptAlein, r_valueAlein, p_valueAlein, std_errAlein = stats.linregress(l_Al * xAl_einohne**2 -xAl_einohne**3/3 , D_Al_einDiff)

plt.plot(l_Al * xAl_einohne**2 -xAl_einohne**3/3 ,D_Al_einDiff , 'bx', label="Messdaten")
x_plot = np.linspace(0, 0.1, 1000)
plt.plot(x_plot,interceptAlein+slopeAlein*x_plot, 'k-', label="Lineare Regression")
plt.yticks( [0 ,1e-03,2e-03, 3e-03, 4e-03, 5e-03],
            [0, 1, 2, 3, 4, 5]
)
plt.legend(loc="best")
plt.xlabel(r'$L-x^2-x^3/3$')
plt.ylabel(r'Durchbiegung D/mm')
plt.grid()
plt.tight_layout
plt.savefig('build/plotAlein.pdf')
plt.close()

#Elastizitätsmodul für Alein

E_Alein = ((m_schraube+m_aufhaeng+m_Alein1)*const.g)/(2*I_Al*slopeAlein)

# Plot für Kupfer doppelseitig belastet ##########################################################################################################################################

slopeCudopl, interceptCudopl, r_valueCudopl, p_valueCudopl, std_errCudopl = stats.linregress(x_dopl(xCu_dopohne[0:7], l_Cu), D_Cu_dopDiff[0:7])

#linke Seite

plt.plot(x_dopl( xCu_dopohne[0:7], l_Cu), D_Cu_dopDiff[0:7] , 'bx', label="Messdaten") # Messpunkte linke Seite
x_plotl = np.linspace(0, 0.3, 1000)
plt.plot(x_plotl,interceptCudopl+slopeCudopl*x_plotl, 'k-', label=r"Lineare Regression $0 \leq x \leq \frac{L}{2} $")
plt.legend(loc="best")
plt.xlabel(r'$3L^2x-4x^3$')
plt.ylabel(r'Durchbiegung D/m')
plt.grid()
plt.tight_layout
plt.savefig('build/plotCudopl.pdf')
plt.close()

# Elastizitätsmodul für Cudop links

E_Cudopl = (const.g*(m_schraube+m_aufhaeng+m_Cudop1+m_Cudop2+m_Cudop3+m_Cudop4))/(48*I_Cu*slopeCudopl)

#rechte Seite

slopeCudopr, interceptCudopr, r_valueCudopr, p_valueCudopr, std_errCudopr = stats.linregress(x_dopr(xCu_dopohne[7:14],l_Cu) , D_Cu_dopDiff[7:14])

plt.plot(x_dopr(xCu_dopohne[7:14], l_Cu) ,D_Cu_dopDiff[7:14] , 'bx', label="Messdaten") # Messpunkte rechte Seite
x_plotr = np.linspace(0, 0.25, 1000)
plt.plot(x_plotr,interceptCudopr+slopeCudopr*x_plotr, 'k-', label=r"Lineare Regression $\frac{L}{2} \leq x \leq L $")

#plt.yticks( [0 ,1e-03,2e-03, 3e-03, 4e-03, 5e-03],
#            [0, 1, 2, 3, 4, 5]
#)
plt.legend(loc="best")
plt.xlabel(r'$4x^3-12Lx^2+9L^2x-L^3$')
plt.ylabel(r'Durchbiegung D/mm')
plt.grid()
plt.tight_layout
plt.savefig('build/plotCudopr.pdf')
plt.close()

# Elastizitätsmodul für Cudop rechts

E_Cudopr = (const.g*(m_schraube+m_aufhaeng+m_Cudop1+m_Cudop2+m_Cudop3+m_Cudop4))/(48*I_Cu*slopeCudopr)

# Plot für Aluminium doppelseitig belastet ##########################################################################################################################################

slopeAldopl, interceptAldopl, r_valueAldopl, p_valueAldopl, std_errAldopl = stats.linregress(x_dopl(xAl_dopohne[0:7], l_Al), D_Al_dopDiff[0:7])

#linke Seite

plt.plot(x_dopl( xAl_dopohne[0:7], l_Al), D_Al_dopDiff[0:7] , 'bx', label="Messdaten") # Messpunkte linke Seite
x_plotl = np.linspace(0, 0.3, 1000)
plt.plot(x_plotl,interceptAldopl+slopeAldopl*x_plotl, 'k-', label=r"Lineare Regression $0 \leq x \leq \frac{L}{2} $")
plt.legend(loc="best")
plt.xlabel(r'$3L^2x-4x^3$')
plt.ylabel(r'Durchbiegung D/m')
plt.grid()
plt.tight_layout
plt.savefig('build/plotAldopl.pdf')
plt.close()

# Elastizitätsmodul für Aldop links

E_Aldopl = (const.g*(m_schraube+m_aufhaeng+m_Aldop1+m_Aldop2+m_Aldop3))/(48*I_Al*slopeCudopl)


#rechte Seite

slopeAldopr, interceptAldopr, r_valueAldopr, p_valueAldopr, std_errAldopr = stats.linregress(x_dopr(xAl_dopohne[7:14],l_Al) , D_Al_dopDiff[7:14])

plt.plot(x_dopr(xAl_dopohne[7:14], l_Al) ,D_Al_dopDiff[7:14] , 'bx', label="Messdaten") # Messpunkte rechte Seite
x_plotr = np.linspace(0, 0.25, 1000)
plt.plot(x_plotr,interceptAldopr+slopeAldopr*x_plotr, 'k-', label=r"Lineare Regression $\frac{L}{2} \leq x \leq L $")

#plt.yticks( [0 ,1e-03,2e-03, 3e-03, 4e-03, 5e-03],
#            [0, 1, 2, 3, 4, 5]
#)
plt.legend(loc="best")
plt.xlabel(r'$4x^3-12Lx^2+9L^2x-L^3$')
plt.ylabel(r'Durchbiegung D/mm')
plt.grid()
plt.tight_layout
plt.savefig('build/plotAldopr.pdf')
plt.close()

# Elastizitätsmodul für Aldop rechts

E_Aldopr = (const.g*(m_schraube+m_aufhaeng+m_Aldop1+m_Aldop2+m_Aldop3))/(48*I_Al*slopeCudopr)

# Tex-Dateien ######################################################################################################

# tex file for Länge des Kupferstabes

with open('build/l_Cu.tex', 'w') as f:
  f.write(make_SI(l_Cu*1e+02,r'\centi\meter', figures=1))

# tex file for Länge des Aluminiumstabes

with open('build/l_Al.tex', 'w') as f:
  f.write(make_SI(l_Al*1e+02,r'\centi\meter', figures=1))

# tex file for Durchmesser des Aluminiumstabes

with open('build/d_Al.tex', 'w') as f:
  f.write(make_SI(np.mean(d_Al)*1e+03,r'\centi\meter', figures=1))

# tex file for Breite des Kupferstabes

with open('build/b_Cu.tex', 'w') as f:
  f.write(make_SI(np.mean(b_Cu)*1e+03,r'\centi\meter', figures=1))

# tex file for Dicke des Kupferstabes

with open('build/d_Cu.tex', 'w') as f:
  f.write(make_SI(np.mean(d_Cu)*1e+03,r'\centi\meter', figures=1))

# tex file for m_aufhaeng

with open('build/m_aufhaeng.tex', 'w') as f:
  f.write(make_SI(m_aufhaeng*1e+03,r'\gram', figures=1))

# tex file for m_schraube
with open('build/m_schraube.tex', 'w') as f:
  f.write(make_SI(m_schraube*1e+03,r'\gram', figures=1))

# tex file for m_Cuein1

with open('build/m_Cuein1.tex', 'w') as f:
  f.write(make_SI(m_Cuein1*1e+03,r'\gram', figures=1))

# tex file for m_Cuein2

with open('build/m_Cuein2.tex', 'w') as f:
  f.write(make_SI(m_Cuein2*1e+03,r'\gram', figures=1))

# tex file for m_Alein1

with open('build/m_Alein1.tex', 'w') as f:
  f.write(make_SI(m_Alein1*1e+03,r'\gram', figures=1))

# tex file for m_Cudop1

with open('build/m_Cudop1.tex', 'w') as f:
  f.write(make_SI(m_Cudop1*1e+03,r'\gram', figures=1))

# tex file for m_Cudop2

with open('build/m_Cudop2.tex', 'w') as f:
  f.write(make_SI(m_Cudop2*1e+03,r'\gram', figures=1))

# tex file for m_Cudop3

with open('build/m_Cudop3.tex', 'w') as f:
  f.write(make_SI(m_Cudop3*1e+03,r'\gram', figures=1))

# tex file for m_Cudop4

with open('build/m_Cudop4.tex', 'w') as f:
  f.write(make_SI(m_Cudop4*1e+03,r'\gram', figures=1))

# tex file for m_Aldop1

with open('build/m_Aldop1.tex', 'w') as f:
  f.write(make_SI(m_Aldop1*1e+03,r'\gram', figures=1))

# tex file for m_Aldop2

with open('build/m_Aldop2.tex', 'w') as f:
  f.write(make_SI(m_Aldop2*1e+03,r'\gram', figures=1))

# tex file for m_Aldop3

with open('build/m_Aldop3.tex', 'w') as f:
  f.write(make_SI(m_Aldop3*1e+03,r'\gram', figures=1))

# tex file for I_Cu

with open('build/I_Cu.tex', 'w') as f:
  f.write(make_SI(I_Cu*1e+10,r'\meter\tothe{4}', exp='1e-10', figures=1))

# tex file for I_Al

with open('build/I_Al.tex', 'w') as f:
  f.write(make_SI(I_Al*1e+10,r'\meter\tothe{4}', exp='1e-10', figures=1))

# tex file for E_Cuein

with open('build/E_Cuein.tex', 'w') as f:
  f.write(make_SI(E_Cuein*1e-09,r'\giga\pascal', exp='1e+09', figures=1))

# tex file for E_Cudopl

with open('build/E_Cudopl.tex', 'w') as f:
  f.write(make_SI(E_Cudopl*1e-09,r'\giga\pascal', exp='1e+09', figures=1))

# tex file for E_Cudopr

with open('build/E_Cudopr.tex', 'w') as f:
  f.write(make_SI(E_Cudopr*1e-09,r'\giga\pascal', exp='1e+09', figures=1))

# tex file for E_Alein

with open('build/E_Alein.tex', 'w') as f:
  f.write(make_SI(E_Alein*1e-09,r'\giga\pascal', exp='1e+09', figures=1))

# tex file for E_Aldopl

with open('build/E_Aldopl.tex', 'w') as f:
  f.write(make_SI(E_Aldopl*1e-09,r'\giga\pascal', exp='1e+09', figures=1))

# tex file for E_Aldopr

with open('build/E_Aldopr.tex', 'w') as f:
  f.write(make_SI(E_Aldopr*1e-09,r'\giga\pascal', exp='1e+09', figures=1))

# Parameter der linearen Regressionen ###############################################################################################################################

# Cuein ---------------------------------------------------------------------------------------------------------------------------------------

# tex file for m from linear regression Cuein

with open('build/m_PlotCuein.tex', 'w') as f:
  f.write(make_SI(slopeCuein,r'', figures=3))

# tex file for b from linear regression Cuein

with open('build/b_PlotCuein.tex', 'w') as f:
  f.write(make_SI(interceptCuein,r'', figures=4))

# Alein ---------------------------------------------------------------------------------------------------------------------------------------

# tex file for m from linear regression Alein

with open('build/m_PlotAlein.tex', 'w') as f:
  f.write(make_SI(slopeAlein,r'', figures=3))

# tex file for b from linear regression Alein

with open('build/b_PlotAlein.tex', 'w') as f:
  f.write(make_SI(interceptAlein,r'', figures=4))

# Cudop ---------------------------------------------------------------------------------------------------------------------------------------

with open('build/m_PlotCudopl.tex', 'w') as f:
  f.write(make_SI(slopeCudopl,r'', figures=3))

# tex file for b from linear regression Aldopl

with open('build/b_PlotCudopl.tex', 'w') as f:
  f.write(make_SI(interceptCudopl,r'', figures=4))

# tex file for m from linear regression Aldopr

with open('build/m_PlotCudopr.tex', 'w') as f:
  f.write(make_SI(slopeCudopr,r'', figures=3))

# tex file for b from linear regression Aldopr

with open('build/b_PlotCudopr.tex', 'w') as f:
  f.write(make_SI(interceptCudopr,r'', figures=4))

# Aldop ---------------------------------------------------------------------------------------------------------------------------------------

# tex file for m from linear regression Aldopl

with open('build/m_PlotAldopl.tex', 'w') as f:
  f.write(make_SI(slopeAldopl,r'', figures=3))

# tex file for b from linear regression Aldopl

with open('build/b_PlotAldopl.tex', 'w') as f:
  f.write(make_SI(interceptAldopl,r'', figures=4))

# tex file for m from linear regression Aldopr

with open('build/m_PlotAldopr.tex', 'w') as f:
  f.write(make_SI(slopeAldopr,r'', figures=3))

# tex file for b from linear regression Aldopr

with open('build/b_PlotAldopr.tex', 'w') as f:
  f.write(make_SI(interceptAldopr,r'', figures=4))

# Tabellen ###############################################################################################################################

# Cu_ein --------------------------------------------------------------------------------------------------------------------------------

xCu_einohne1, xCu_einohne2 = np.array_split(xCu_einohne,2)
DCu_einohne1, DCu_einohne2 = np.array_split(DCu_einohne,2)
DCu_einmit1, DCu_einmit2 = np.array_split(DCu_einmit,2)
D_Cu_einDiff1, D_Cu_einDiff2 = np.array_split(D_Cu_einDiff,2)

table_header = r'''
  \begin{tabular}{c c c c c c c c}
    \toprule
     & \multicolumn{1}{c}{$D_0$} & \multicolumn{1}{c}{$D_m$} & \multicolumn{1}{c}{Differenz} & & \multicolumn{1}{c}{$D_0$} & \multicolumn{1}{c}{$D_m$} & \multicolumn{1}{c}{Differenz}\\
    \cmidrule(lr{0,5em}){1-4} \cmidrule(lr{0,5em}){5-8}
    {$x \:/\: \si{\centi\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} &
    {$x \:/\: \si{\centi\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} \\

    \cmidrule(lr{0,5em}){1-4} \cmidrule(lr{0,5em}){5-8}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.2f} & {2:1.2f} & {3:1.2f} & {4:1.1f} & {5:1.2f} & {6:1.2f} & {7:1.2f} \\'

with open('build/Cu_ein.tex', 'w') as g:
    g.write(table_header)
    for row in zip(xCu_einohne1*1e+02, DCu_einohne1, DCu_einmit1, D_Cu_einDiff1*1e+03, xCu_einohne2*1e+02, DCu_einohne2, DCu_einmit2, D_Cu_einDiff2*1e+03):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

# Cu_dop --------------------------------------------------------------------------------------------------------------------------------

xCu_dopohne1, xCu_dopohne2 = np.array_split(xCu_dopohne,2)
DCu_dopohne1, DCu_dopohne2 = np.array_split(DCu_dopohne,2)
DCu_dopmit1, DCu_dopmit2 = np.array_split(DCu_dopmit,2)
D_Cu_dopDiff1, D_Cu_dopDiff2 = np.array_split(D_Cu_dopDiff,2)

table_header = r'''
  \begin{tabular}{c c c c c c c c}
    \toprule
     & \multicolumn{1}{c}{$D_0$} & \multicolumn{1}{c}{$D_m$} & \multicolumn{1}{c}{Differenz} & & \multicolumn{1}{c}{$D_0$} & \multicolumn{1}{c}{$D_m$} & \multicolumn{1}{c}{Differenz}\\
    \cmidrule(lr{0,5em}){1-4} \cmidrule(lr{0,5em}){5-8}
    {$x \:/\: \si{\centi\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} &
    {$x \:/\: \si{\centi\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} \\

    \cmidrule(lr{0,5em}){1-4} \cmidrule(lr{0,5em}){5-8}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.2f} & {2:1.2f} & {3:1.2f} & {4:1.1f} & {5:1.2f} & {6:1.2f} & {7:1.2f} \\'

with open('build/Cu_dop.tex', 'w') as g:
    g.write(table_header)
    for row in zip(xCu_dopohne1*1e+02, DCu_dopohne1, DCu_dopmit1, D_Cu_dopDiff1*1e+03, xCu_dopohne2*1e+02, DCu_dopohne2, DCu_dopmit2, D_Cu_dopDiff2*1e+03):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

# Al_ein --------------------------------------------------------------------------------------------------------------------------------

xAl_einohne1, xAl_einohne2 = np.array_split(xAl_einohne,2)
DAl_einohne1, DAl_einohne2 = np.array_split(DAl_einohne,2)
DAl_einmit1, DAl_einmit2 = np.array_split(DAl_einmit,2)
D_Al_einDiff1, D_Al_einDiff2 = np.array_split(D_Al_einDiff,2)

table_header = r'''
  \begin{tabular}{c c c c c c c c}
    \toprule
     & \multicolumn{1}{c}{$D_0$} & \multicolumn{1}{c}{$D_m$} & \multicolumn{1}{c}{Differenz} & & \multicolumn{1}{c}{$D_0$} & \multicolumn{1}{c}{$D_m$} & \multicolumn{1}{c}{Differenz}\\
    \cmidrule(lr{0,5em}){1-4} \cmidrule(lr{0,5em}){5-8}
    {$x \:/\: \si{\centi\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} &
    {$x \:/\: \si{\centi\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} \\

    \cmidrule(lr{0,5em}){1-4} \cmidrule(lr{0,5em}){5-8}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.2f} & {2:1.2f} & {3:1.2f} & {4:1.1f} & {5:1.2f} & {6:1.2f} & {7:1.2f} \\'

with open('build/Al_ein.tex', 'w') as g:
    g.write(table_header)
    for row in zip(xAl_einohne1*1e+02, DAl_einohne1, DAl_einmit1, D_Al_einDiff1*1e+03, xAl_einohne2*1e+02, DAl_einohne2, DAl_einmit2, D_Al_einDiff2*1e+03):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

# Al_dop --------------------------------------------------------------------------------------------------------------------------------

xAl_dopohne1, xAl_dopohne2 = np.array_split(xAl_dopohne,2)
DAl_dopohne1, DAl_dopohne2 = np.array_split(DAl_dopohne,2)
DAl_dopmit1, DAl_dopmit2 = np.array_split(DAl_dopmit,2)
D_Al_dopDiff1, D_Al_dopDiff2 = np.array_split(D_Al_dopDiff,2)

table_header = r'''
  \begin{tabular}{c c c c c c c c}
    \toprule
     & \multicolumn{1}{c}{$D_0$} & \multicolumn{1}{c}{$D_m$} & \multicolumn{1}{c}{Differenz} & & \multicolumn{1}{c}{$D_0$} & \multicolumn{1}{c}{$D_m$} & \multicolumn{1}{c}{Differenz}\\
    \cmidrule(lr{0,5em}){1-4} \cmidrule(lr{0,5em}){5-8}
    {$x \:/\: \si{\centi\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} &
    {$x \:/\: \si{\centi\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} & {$D(x) \:/\: \si{\micro\meter}$} \\

    \cmidrule(lr{0,5em}){1-4} \cmidrule(lr{0,5em}){5-8}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.2f} & {2:1.2f} & {3:1.2f} & {4:1.1f} & {5:1.2f} & {6:1.2f} & {7:1.2f} \\'

with open('build/Al_dop.tex', 'w') as g:
    g.write(table_header)
    for row in zip(xAl_dopohne1*1e+02, DAl_dopohne1, DAl_dopmit1, D_Al_dopDiff1*1e+03, xAl_dopohne2*1e+02, DAl_dopohne2, DAl_dopmit2, D_Al_dopDiff2*1e+03):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)


# Testprints ###############################################################################################################################

print(interceptAlein)
print(I_Cu)
print(I_Al)
print(E_Alein)
print(m_Al_stange/(np.pi*(np.mean(d_Al/2)**2)*l_Al))
print(E_Alein)
print(E_Aldopl)
print(E_Aldopr)