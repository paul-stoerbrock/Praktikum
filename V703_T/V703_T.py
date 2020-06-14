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


def rel_err(mess, lit):
    return abs((unp.noms(mess)-unp.noms(lit))/unp.noms(lit))*100



# Messdaten #################################################################

U_kenn, N = np.genfromtxt('Kennlinie.dat', unpack=True)

U_zaehl, I = np.genfromtxt('Zaehlrohrstrom.dat', unpack=True)


N_err = np.sqrt(N)

N_miterr = unp.uarray(N, N_err)

I_miterr = unp.uarray(I*1e-06, 0.05e-06)

# Plots ######################################################################

# Plot zur Kennlinie des Geiger-Müller Zählrohrs

par,cov = np.polyfit(U_kenn[3:36], N[3:36], deg=1, cov=True)
err = np.sqrt(np.diag(cov))

plt.errorbar(U_kenn, N, yerr=stds(N_miterr) ,fmt='kx', label='Messwerte mit Fehler')
x_plot = np.linspace(350, 680, 10000)
plt.plot(x_plot,par[0]*x_plot+par[1] , 'r-', label='Fitkurve')
plt.legend(loc="best")
plt.xlabel(r'Spannung $U \:/\:V$')
plt.ylabel(r'Intensität $I\:/\:Imp/60s$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_kenn.pdf')
plt.close()

par_kenn=unp.uarray(par/60,err/60)

m_winkel = unp.arctan(par_kenn[0]/100)
m_percent = (m_winkel*180)/np.pi

# Plot zu Z #################################################################

N_Z = np.array([9837, 9995, 10264, 10151, 10184, 10253, 10493, 11547])
N_Z_miterr= unp.uarray(N_Z, np.sqrt(N_Z))
Z = I_miterr/(const.e*N_Z_miterr)

par, cov= np.polyfit(noms(I_miterr) ,noms(Z) ,deg=1 ,cov=True)
err = np.sqrt(np.diag(cov))

plt.errorbar(noms(I_miterr) ,noms(Z) , xerr=(stds(I_miterr)) ,yerr=stds(Z) , fmt='kx', label='Messwerte mit Fehler')
x_plot = np.linspace(0.2*1e-06, 1.8*1e-06, 10000)
plt.plot(x_plot,par[0]*x_plot+par[1] , 'r-', label='Fitkurve')
plt.xticks([0.2*1e-06, 0.4*1e-06, 0.6*1e-06, 0.8*1e-06, 1.0*1e-06, 1.2*1e-06, 1.4*1e-06, 1.6*1e-06, 1.8*1e-06],
[0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8])
plt.legend(loc="best")
plt.ylabel(r'freigesetzte Ladung pro einfallendem Teilchen Z ')
plt.xlabel(r'Stromstärke $I\:/\:\mu A $')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_Z.pdf')
plt.close()

par_Z=unp.uarray(par, err)

#Berechnungen Zwei Quellen Methode

N1 = ufloat(96041, np.sqrt(96401))
N2 = ufloat(76518,np.sqrt(76518))
N12= ufloat(158479, np.sqrt(158479))

T = (N1+N2-N12)/(2*N1*N2) *120



# Tex Files ######################################################################

with open('build/a_kenn.tex', 'w') as f: 
  f.write(make_SI(par_kenn[0]*1e02 ,r'\second\tothe{-1}\volt\tothe{-1} ',exp='e-02' ,figures=2))

with open('build/b_kenn.tex', 'w') as f: 
  f.write(make_SI(par_kenn[1] ,r'\second\tothe{-1} ' ,figures=2))


# Tex Files of Totzeit

with open('build/N_1.tex', 'w') as f: 
  f.write(make_SI(N1 ,r'Imp/120\second' ,figures=2))

with open('build/N_2.tex', 'w') as f: 
  f.write(make_SI(N2 ,r'Imp/120\second ' ,figures=2))

with open('build/N_12.tex', 'w') as f: 
  f.write(make_SI(N12 ,r'Imp/120\second ' ,figures=2))

with open('build/T.tex', 'w') as f: 
  f.write(make_SI(T*1e06 ,r'\micro\second ' ,figures=2))



with open('build/a_Z.tex', 'w') as f: 
  f.write(make_SI(par_Z[0]*1e-14 ,r'\ampere\tothe{-1} ',exp='e14' ,figures=2))

with open('build/b_Z.tex', 'w') as f: 
  f.write(make_SI(par_Z[1]*1e-7 ,r' ',exp='e07' ,figures=2))

with open('build/m.tex', 'w') as f: 
  f.write(make_SI(m_percent*100 ,r'\percent\per{100}\volt ' ,figures=2))




# Tabellen ####################################################################

N1, N2, N3= np.array_split(N,3)
U1, U2, U3= np.array_split(U_kenn,3)

table_header = r'''
  \begin{longtable}[H]{S[table-format=3.0] S[table-format=5.0]@{${}\pm{}$} S[table-format=3.0] S[table-format=3.0] S[table-format=5.0]@{${}\pm{}$} S[table-format=3.0] S[table-format=3.0] S[table-format=5.0]@{${}\pm{}$} S[table-format=3.0]}
    \caption{Messwerte für die Charakteristik des Zählrohrs.
    Die Messung ist in Abständen von \SI{10}{\volt} durchgeführt worden und
    der Fehler bei N kommt von der Poisson-Verteilung $ \Delta N = \sqrt{N}$.
    }\\
    \toprule
    \multicolumn{1}{c}{ $U\:/\:V$ } & \multicolumn{2}{c}{$N\:/\:Imp/60s$ }   &
    \multicolumn{1}{c}{ $U\:/\:V$ } & \multicolumn{2}{c}{  $N\:/\:Imp/60s$ }   &
    \multicolumn{1}{c}{ $U\:/\:V$ } & \multicolumn{2}{c}{  $N\:/\:Imp/60s$ } \\
    \cmidrule(lr{0,5em}){1-3} \cmidrule(lr{0,5em}){4-6}  \cmidrule(lr{0,5em}){7-9}

'''
table_footer = r'''    \bottomrule
  \label{tab:1}
  \end{longtable}
'''
row_template = r'     {0:1.0f} & {1:1.0f} & {2:1.0f} & {3:1.0f} & {4:1.0f} & {5:1.0f}& {6:1.0f} & {7:1.0f} & {8:1.0f}   \\'

with open('build/kenn_table.tex', 'w') as g:
    g.write(table_header)
    for row in zip( U1, N1, np.sqrt(N1), U2, N2, np.sqrt(N2), U3, N3, np.sqrt(N3)):
        g.write(row_template.format(*row).replace('nan',' ').replace('.',',') )
        g.write('\n')
    g.write(table_footer)


print(Z)

table_header = r'''
  \begin{tabular}[H]{S[table-format=3.0] S[table-format=5.0] @{${}\pm{}$} S[table-format=3.0] S[table-format=1.1]@{${}\pm{}$}S[table-format=0.2]  S[table-format=1.2]@{${}\pm{}$}S[table-format=0.2]}
    \toprule
    \multicolumn{1}{c}{ $U\:/\:V$ } & \multicolumn{2}{c}{$N\:/\:Imp/60s$ }   &
    \multicolumn{2}{c}{ $I\:/\:\mu A$ } & \multicolumn{2}{c}{Z /$10^8$} \\
    \cmidrule(lr{0,5em}){1-7}

'''
table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.0f} & {1:1.0f} & {2:1.0f} & {3:1.1f} & {4:1.2f} & {5:1.2f} & {6:1.2f} \\'

with open('build/Z_table.tex', 'w') as g:
    g.write(table_header)
    for row in zip(U_zaehl ,noms(N_Z_miterr), stds(N_Z_miterr), noms(I_miterr)*1e06, stds(I_miterr)*1e06,noms(Z)*1e-08, stds(Z)*1e-08):
        g.write(row_template.format(*row).replace('nan',' ').replace('.',',') )
        g.write('\n')
    g.write(table_footer)