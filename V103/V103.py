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

# Plot für Kupfer einseitig belastet

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

# Für Cuein

E_Cuein = ((m_aufhaeng+m_schraube+m_Cuein1+m_Cuein2)*9.81)/(2*I_Cu*slopeCuein)

# Plot für Aluminium einseitig belastet ####################################################################################################################

# Für Alein

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

# Für Alein

E_Alein = ((m_schraube+m_aufhaeng+m_Alein1)*9.81)/(2*I_Al*slopeAlein)

# Plot für Kupfer doppelseitig belastet ##########################################################################################################################################

slopeCudopl, interceptCudopl, r_valueCudopl, p_valueCudopl, std_errCudopl = stats.linregress(3*l_Cu**2*xCu_dopohne[0:7]-4*xCu_dopohne[0:7]**3 , D_Cu_dopDiff[0:7])

#linke Seite

plt.plot(3*l_Cu**2*xCu_dopohne[0:7]-4*xCu_dopohne[0:7]**3 ,D_Cu_einDiff[0:7] , 'bx', label="Messdaten") # Messpunkte linke Seite
x_plotl = np.linspace(0, l_Cu**3, 1000)
plt.plot(x_plotl,interceptCudopl+slopeCudopl*x_plotl, 'k-', label=r"Lineare Regression $0 \leq x \leq \frac{L}{2} $")
plt.legend(loc="best")
plt.xlabel(r'$3L^2x-4x^3$')
plt.ylabel(r'Durchbiegung D/mm')
plt.grid()
plt.tight_layout
plt.savefig('build/plotCudopl.pdf')
plt.close()

#rechte Seite

slopeCudopr, interceptCudopr, r_valueCudopr, p_valueCudopr, std_errCudopr = stats.linregress(4*xCu_dopohne[8:14]**3 - 12 * l_Cu * xCu_dopohne[8:14]**2 +9* l_Cu**2* xCu_dopohne[8:14] - l_Cu**3 , D_Cu_dopDiff[8:14])

plt.plot(4*xCu_dopohne[8:14]**3 - 12 * l_Cu * xCu_dopohne[8:14]**2 + 9 * l_Cu**2 * xCu_dopohne[8:14] - l_Cu**3 ,D_Cu_einDiff[8:14] , 'bx', label="Messdaten") # Messpunkte rechte Seite
x_plotr = np.linspace(l_Cu/2, 0 , 1000)
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
row_template = r'     {0:1.2f} & {1:1.2f} & {2:1.2f} & {3:1.2f} & {4:1.2f} & {5:1.2f} & {6:1.2f} & {7:1.2f} \\'

with open('build/Cu_ein.tex', 'w') as g:
    g.write(table_header)
    for row in zip(xCu_einohne1, DCu_einohne1, DCu_einmit1, D_Cu_einDiff1*1e+03, xCu_einohne2, DCu_einohne2, DCu_einmit2, D_Cu_einDiff2*1e+03):
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
row_template = r'     {0:1.2f} & {1:1.2f} & {2:1.2f} & {3:1.2f} & {4:1.2f} & {5:1.2f} & {6:1.2f} & {7:1.2f} \\'

with open('build/Cu_dop.tex', 'w') as g:
    g.write(table_header)
    for row in zip(xCu_dopohne1, DCu_dopohne1, DCu_dopmit1, D_Cu_dopDiff1*1e+03, xCu_dopohne2, DCu_dopohne2, DCu_dopmit2, D_Cu_dopDiff2*1e+03):
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
row_template = r'     {0:1.2f} & {1:1.2f} & {2:1.2f} & {3:1.2f} & {4:1.2f} & {5:1.2f} & {6:1.2f} & {7:1.2f} \\'

with open('build/Al_ein.tex', 'w') as g:
    g.write(table_header)
    for row in zip(xAl_einohne1, DAl_einohne1, DAl_einmit1, D_Al_einDiff1*1e+03, xAl_einohne2, DAl_einohne2, DAl_einmit2, D_Al_einDiff2*1e+03):
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
row_template = r'     {0:1.2f} & {1:1.2f} & {2:1.2f} & {3:1.2f} & {4:1.2f} & {5:1.2f} & {6:1.2f} & {7:1.2f} \\'

with open('build/Al_dop.tex', 'w') as g:
    g.write(table_header)
    for row in zip(xAl_dopohne1, DAl_dopohne1, DAl_dopmit1, D_Al_dopDiff1*1e+03, xAl_dopohne2, DAl_dopohne2, DAl_dopmit2, D_Al_dopDiff2*1e+03):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)


# Testprints ###############################################################################################################################

print(I_Cu)
print(E_Alein)
print(m_Al_stange/(np.pi*(np.mean(d_Al/2)**2)*l_Al))
print(slopeCudopr*l_Cu**3)
print(D_Al_dopDiff)