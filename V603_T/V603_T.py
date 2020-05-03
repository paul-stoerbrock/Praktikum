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
        x = '{0:.{1:}f}'.format(num , figures)

    return r'\SI{{{}{}}}{{{}}}'.format(x, exp, unit)


# Messwerte ##############################################################################

# Emission von Kupfer
theta_deg_Cu, N_Cu = np.genfromtxt('EmissionCu.dat', unpack=True)

# Zählrate der Röntgenstrahlung mit Al-Absorber
theta_deg_Al, R_Al = np.genfromtxt('ComptonAl.txt', unpack=True)

# Zählrate der Röntgenstrahlung ohne Absorber
theta_deg_Ohne, R_Ohne = np.genfromtxt('ComptonOhne.txt', unpack=True)

# Umrechnung theta zu lambda

theta_rad_Al = np.deg2rad(theta_deg_Al)


lambda1=2 *201.4*1e-12*np.sin(theta_rad_Al)





# Plots ###########################################################################################

# Plot des Kupfer Emissionsspektrums

plt.plot(theta_deg_Cu, N_Cu,'kx', label='Messwerte')
x_plot = np.linspace(0, 5, 1000)
plt.legend(loc="best")
plt.annotate('$K_{\\beta}$', xy=(20, 1700), size = 15

)
plt.annotate('$K_{\\alpha} $', xy=(23, 5050), size = 15

)
plt.annotate('Bremsspektrum', xy=(10, 700), size = 20

)
plt.xlabel(r'Winkel $\theta \:/\:°$')
plt.ylabel(r'Intensität $I\:/\:Imp/s$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_Cu.pdf')
plt.close()



 #Plot des Absorptionsspektrum ohne Absorber




#Totzeitkorrektur


I_Al= R_Al/(1-90*1e-06*R_Al)


I_Ohne= R_Ohne/(1-90*1e-06*R_Ohne)


T = I_Al/I_Ohne

par, covm = np.polyfit(lambda1, unp.nominal_values(T), deg=1, cov=True)

err = np.sqrt(np.diag(covm))



plt.plot(lambda1*1e12, unp.nominal_values(T),'kx', label='Messwerte')
x_plot = np.linspace(4.7*1e-11, 7*1e-11, 1000)
plt.plot(x_plot*1e12,par[0]*x_plot+par[1] ,'r-', label="Emissionsspektrum Cu")
plt.legend(loc="best")
plt.xlabel(r'Wellenlänge $\lambda \:/\: pm$')
plt.ylabel(r'Transmission T')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_T.pdf')
plt.close()

par=unp.uarray(par, err)

I_0 = 2731
I_1 = 1180
I_2 = 1024

print((I_1/I_0-par[1])/par[0] )
print((I_2/I_0-par[1])/par[0])

# Tex Files ##################################################################################################################

#tex file of a_Al  ###########################################################################

with open('build/a_Al.tex', 'w') as f: 
  f.write(make_SI(par[0] ,r'\meter\tothe{-1}' ,figures=1))

#tex file of b_Al  ###########################################################################

with open('build/b_Al.tex', 'w') as f: 
  f.write(make_SI(par[1] ,r'' ,figures=1))




# Tabellen ############################################################################################################################

# Tabelle des Transmissionsspektrum von Aluminium

table_header = r'''
  \begin{tabular}{c c c }
    \toprule
    \multicolumn{1}{c}{Winkel $\theta \:/\:° $} & \multicolumn{1}{c}{Zählrate $R_{Al}\:/\:Imp/s $} & \multicolumn{1}{c}{Zählrate $R_{Ohne}\:/\:Imp/s $} \\
    \cmidrule(lr){1-3} 
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.1f} & {2:1.1f}  \\'


with open('build/table_Al.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_deg_Al, R_Al, R_Ohne ):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)

# Tabelle für das Emissionspektrum von Kupfer

theta_deg_Cu_1, theta_deg_Cu_2, theta_deg_Cu_3 = np.array_split(theta_deg_Cu, 3)
N_Cu_1, N_Cu_2, N_Cu_3 = np.array_split(N_Cu, 3)


table_header = r'''
  \begin{tabular}{c c c c c c}
    \toprule
    \multicolumn{1}{c}{Winkel $\theta \:/\:° $ } & \multicolumn{1}{c}{Zählrate $N\:/\:Imp/s $ } & \multicolumn{1}{c}{Winkel $\theta \:/\:°$ } & \multicolumn{1}{c}{Zählrate $N\:/\:Imp/s $ }& \multicolumn{1}{c}{Winkel $\theta \:/\:° $}& \multicolumn{1}{c}{Zählrate $N\:/\:Imp/s $  }\\
    \cmidrule(lr){1-2} \cmidrule(lr{0,5em}){3-4} \cmidrule(lr{0,5em}){5-6}
'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.0f} & {2:1.1f} & {3:1.0f} & {4:1.1f} & {5:1.0f} \\'


with open('build/table_Cu.tex', 'w') as g:
    g.write(table_header)
    for row in zip(theta_deg_Cu_1, N_Cu_1, theta_deg_Cu_2, N_Cu_2, theta_deg_Cu_3, N_Cu_3):
        g.write(row_template.format(*row))
        g.write('\n')
    g.write(table_footer)