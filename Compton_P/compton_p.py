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





Theta_Cu, N_Cu = np.genfromtxt('EmissionCu.dat', unpack = True)

Alpha_Al, N_Al = np.genfromtxt('ComptonAl.txt', unpack = True)

Alpha_ohne, N_ohne = np.genfromtxt('ComptonOhne.txt', unpack = True)

d_LiF = 201.4*1e-12

tau = 90*1e-06

Lambda = 2*d_LiF*np.sin(np.deg2rad(Alpha_Al))

I_Al = N_Al/(1-tau*N_Al)
I_ohne = N_ohne/(1-tau*N_ohne)

T = I_Al/I_ohne

par, covm = np.polyfit(Theta_Cu, N_Cu, deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(Theta_Cu, N_Cu, 'kx', label='Messwerte')
x_plot = np.linspace(8, 30, 1000)

plt.xlabel(r'Winkel $\alpha/\degree$')
plt.ylabel(r'Impuls/Sekunde')
ymax = np.max(N_Cu)
xmax = Theta_Cu[145]
plt.annotate('Global max', xy=(xmax, ymax), xytext=(22.75, ymax), horizontalalignment='left')
ylmax = N_Cu[122]
xlmax = Theta_Cu[122]
plt.annotate('Local max', xy=(xlmax, ylmax), xytext=(20, ylmax), horizontalalignment='right')
plt.axvline(x=xmax, color='r', linestyle= ':', label='$K_{\\alpha}$-Linie')
plt.axvline(x=xlmax, color='b', linestyle= ':', label='$K_{\\beta}$-Linie')
plt.annotate('$K_{\\alpha}$', xy=(xmax, 2500), xycoords='data', xytext=(xmax+1, 2500), textcoords='data', arrowprops=dict(arrowstyle='->', facecolor='grey'), horizontalalignment='left')
plt.annotate('$K_{\\beta}$', xy=(xlmax, 2500), xycoords='data', xytext=(xlmax+1, 2500), textcoords='data', arrowprops=dict(arrowstyle='->', facecolor='grey'), horizontalalignment='left')
plt.annotate('Bremsstrahlung', xy=(xmax/2, 1100), xytext=(xmax/2, 1100), verticalalignment='top', horizontalalignment='left')
plt.legend(loc="best")
plt.grid()
plt.tight_layout
plt.savefig('build/plotCu.pdf')
plt.close()

xlmax = 20.2  #degrees
xmax = 22.5    #degrees

L_a = 2*d_LiF*np.sin(np.deg2rad(xmax))
L_b = 2*d_LiF*np.sin(np.deg2rad(xlmax))

E_a = (const.h*const.c)/(L_a*const.e)
E_b = (const.h*const.c)/(L_b*const.e)

print(L_a)
print(L_b)
print(E_a)
print(E_b)


par, covm = np.polyfit(Lambda, T, deg=1, cov=True)
err = np.sqrt(np.diag(covm))

plt.plot(Lambda, T, 'kx', label='Messwerte')
x_plot = np.linspace(4.7*1e-11, 7*1e-11, 1000)
plt.plot(x_plot ,par[0]*x_plot+par[1], 'b-', label="Lineare Regression")
plt.legend(loc="best")
plt.xlabel(r'Transmission')
plt.ylabel(r'Wellenlänge $\lambda / m$')
plt.grid()
plt.tight_layout
plt.savefig('build/plotAl.pdf')
plt.close()

m_Al = par[0]
b_Al =par[1]

I_0 = 2731
I_1 = 1180
I_2 = 1024

T_1 = I_1/I_0
T_2 = I_2/I_0

L_1 = (T_1-par[1])/par[0]
L_2 = (T_2-par[1])/par[0]

L_C = L_2-L_1

print(L_1)
print(L_2)
print(L_C)

# tex file for xmax

with open('build/xmax.tex', 'w') as f:
  f.write(make_SI(xmax,r'°', figures=2))

# tex file for xlmax

with open('build/xlmax.tex', 'w') as f:
  f.write(make_SI(xlmax,r'°', figures=2))

# tex file for L_a

with open('build/L_a.tex', 'w') as f:
  f.write(make_SI(L_a*1e+10,r'\angstrom', figures=2))

# tex file for L_b

with open('build/L_b.tex', 'w') as f:
  f.write(make_SI(L_b*1e+10,r'\angstrom', figures=2))

# tex file for E_a

with open('build/E_a.tex', 'w') as f:
  f.write(make_SI(E_a,r'\electronvolt', figures=2))

# tex file for E_b

with open('build/E_b.tex', 'w') as f:
  f.write(make_SI(E_b,r'\electronvolt', figures=2))

# tex file for m_Al

with open('build/m_Al.tex', 'w') as f:
  f.write(make_SI(par[0]*1e-10,r'\meter', exp='1e+10', figures=1))

# tex file for b_Al

with open('build/b_Al.tex', 'w') as f:
  f.write(make_SI(par[1],r'', figures=2))

# tex file for I_0

with open('build/I_0.tex', 'w') as f:
  f.write(make_SI(I_0,r'Impuls\per\second', figures=0))

# tex file for I_1

with open('build/I_1.tex', 'w') as f:
  f.write(make_SI(I_1,r'Impuls\per\second', figures=0))

# tex file for I_2

with open('build/I_2.tex', 'w') as f:
  f.write(make_SI(I_2,r'\impuls\per\second', figures=0))

# tex file for T_1

with open('build/T_1.tex', 'w') as f:
  f.write(make_SI(T_1,r'', figures=2))

# tex file for T_2

with open('build/T_2.tex', 'w') as f:
  f.write(make_SI(T_2,r'', figures=2))

# tex file for L_1

with open('build/L_1.tex', 'w') as f:
  f.write(make_SI(L_1,r'\meter', figures=2))

# tex file for L_2

with open('build/L_2.tex', 'w') as f:
  f.write(make_SI(L_2,r'\meter', figures=2))

# tex file for L_C

with open('build/L_C.tex', 'w') as f:
  f.write(make_SI(L_C*1e+12,r'\meter', exp='1e-12', figures=2))



Theta_Cu1, Theta_Cu2, Theta_Cu3 = np.array_split(Theta_Cu[0:99],3)
N_Cu1, N_Cu2, N_Cu3 = np.array_split(N_Cu[0:99],3)

table_header = r'''
  \begin{tabular}{c c c c c c}
    \toprule
    \multicolumn{1}{c}{$\text{Winkel} \; \alpha\:/\si{\degree}$} & \multicolumn{1}{c}{$\text{Impuls/Sekunde}$} & 
    \multicolumn{1}{c}{$\text{Winkel} \; \alpha\:/\si{\degree}$} & \multicolumn{1}{c}{$\text{Impuls/Sekunde}$} &  
    \multicolumn{1}{c}{$\text{Winkel} \; \alpha\:/\si{\degree}$} & \multicolumn{1}{c}{$\text{Impuls/Sekunde}$}\\

    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4} \cmidrule(lr{0,5em}){5-6}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.2f} & {2:1.2f} & {3:1.2f} & {4:1.2f} & {5:1.2f}\\'

with open('build/table_Cu1.tex', 'w') as g:
    g.write(table_header)
    for row in zip(Theta_Cu1, N_Cu1, Theta_Cu2, N_Cu2, Theta_Cu3, N_Cu3):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

Theta_Cu1, Theta_Cu2, Theta_Cu3 = np.array_split(Theta_Cu[100:170],3)
N_Cu1, N_Cu2, N_Cu3 = np.array_split(N_Cu[100:170],3)

table_header = r'''
  \begin{tabular}{c c c c c c}
    \toprule
    \multicolumn{1}{c}{$\text{Winkel} \; \alpha\:/\si{\degree}$} & \multicolumn{1}{c}{$\text{Impuls/Sekunde}$} & 
    \multicolumn{1}{c}{$\text{Winkel} \; \alpha\:/\si{\degree}$} & \multicolumn{1}{c}{$\text{Impuls/Sekunde}$} &  
    \multicolumn{1}{c}{$\text{Winkel} \; \alpha\:/\si{\degree}$} & \multicolumn{1}{c}{$\text{Impuls/Sekunde}$}\\

    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4} \cmidrule(lr{0,5em}){5-6}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.2f} & {2:1.2f} & {3:1.2f} & {4:1.2f} & {5:1.2f}\\'

with open('build/table_Cu2.tex', 'w') as g:
    g.write(table_header)
    for row in zip(Theta_Cu1, N_Cu1, Theta_Cu2, N_Cu2, Theta_Cu3, N_Cu3):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)


  
Alpha_Al1, Alpha_Al2 = np.array_split(Alpha_Al,2)
N_Al1, N_Al2 = np.array_split(N_Al,2)
N_ohne1, N_ohne2 = np.array_split(N_ohne,2)


table_header = r'''
  \begin{tabular}{c c c c c c}
    \toprule
    \multicolumn{1}{c}{$\text{Winkel} \; \beta\:/\si{\degree}$} & \multicolumn{1}{c}{$\text{$N_{AL}$/Sekunde}$} & \multicolumn{1}{c}{$\text{$N_0$/Sekunde}$} &
    \multicolumn{1}{c}{$\text{Winkel} \; \beta\:/\si{\degree}$} & \multicolumn{1}{c}{$\text{$N_{Al}$/Sekunde}$} & \multicolumn{1}{c}{$\text{$N_0$/Sekunde}$}\\
    \cmidrule(lr{0,5em}){1-3} \cmidrule(lr{0,5em}){4-6}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.2f} & {2:1.2f} & {3:1.2f} & {4:1.2f} & {5:1.2f}\\'

with open('build/table_Al.tex', 'w') as g:
    g.write(table_header)
    for row in zip(Alpha_Al1, N_Al1, N_ohne1, Alpha_Al2, N_Al2, N_ohne2):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)
