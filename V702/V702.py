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
import sympy

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
    return abs((unp.nominal_values(mess)-unp.nominal_values(lit))/unp.nominal_values(lit))*100

def error(f, err_vars=None):
    from sympy import Symbol, latex
    s = 0
    latex_names = dict()
    
    if err_vars == None:
        err_vars = f.free_symbols
        
    for v in err_vars:
        err = Symbol('latex_std_' + v.name)
        s += f.diff(v)**2 * err**2
        latex_names[err] = '\\sigma_{' + latex(v) + '}'
        
    return latex(sympy.sqrt(s), symbol_names=latex_names)


def f(x, A, B, C, D):
    return A*(np.exp(B*x)+np.exp(C*x))+D



# Messwerte #########################################################

t_V, N_V = np.genfromtxt('Vanadium.dat' , unpack=True)
t_Rh, N_Rh = np.genfromtxt('Rhodium.dat', unpack=True)

N_U_V= np.array([129,129,129,129,129,129,129,129,129,129,143,143,143,143,143,143,143,143,143,143,144,144,144,144,144,144,144,144,144,144,136,136,136,136,136,136,136,136,136,136,139,139,139,139])
N_U_V_err = unp.uarray(N_U_V, np.sqrt(N_U_V))

N_V_err = unp.uarray(N_V, np.sqrt(N_V))

N_V_ohne_U = N_V_err-N_U_V_err/10



N_U_Rh =np.array([129,129,129,129,129,129,129,129,129,129,129,129,129,129,129,129,129,129,129,129,143,143,143,143,143,143,143,143,143,143,143,143,143,143,143,143,143,143,143,143,136,136,136,136])
N_U_Rh_err =unp.uarray(N_U_Rh,np.sqrt(N_U_Rh))

N_Rh_err = unp.uarray(N_Rh, np.sqrt(N_Rh))

N_Rh_ohne_U = N_Rh_err-N_U_Rh_err/20




# Plots #####################################################################

# Plot von Vanadium %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par, cov= np.polyfit(t_V, np.log(noms(N_V_ohne_U)), deg=1, cov=True)
err= np.sqrt(np.diag(cov))

plt.errorbar(t_V, noms(N_V_ohne_U),xerr=stds(N_V_ohne_U) ,fmt='ko', label='Messwerte')
x_plot = np.linspace(0, 1400, 10000)
plt.plot(x_plot, np.exp(x_plot*par[0]+par[1]), 'r-', label='Ausgleichsgerade')
plt.yscale('log')
plt.legend(loc="best")
plt.xlabel(r' $t \:/\:s$')
plt.ylabel(r' $I\:/\:Imp/30s$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_V.pdf')
plt.close()

par_V = unp.uarray(par, err)

# Plot von Rhodium %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par, cov = np.polyfit( t_Rh[18:], np.log(noms(N_Rh_ohne_U[18:])), deg=1, cov=True)
err1 = np.sqrt(np.diag(cov))


par3, cov3 = np.polyfit(t_Rh[:17], np.log(noms(N_Rh_ohne_U[:17])- np.exp(par[1]+par[0]*t_Rh[:17])  ), deg=1, cov=True)
err3 = np.sqrt(np.diag(cov3))

plt.errorbar(t_Rh, noms(N_Rh_ohne_U),xerr=stds(N_Rh_ohne_U) ,fmt='kx', label='Messwerte')
x_plot = np.linspace(0, 660, 10000)
x_plot2 = np.linspace(0, 300, 1000)
plt.plot(x_plot, np.exp(x_plot*par[0]+par[1]) , 'r-', label='Gerade des langlebigen Zerfalls')
plt.plot(x_plot,np.exp(x_plot*par[0]+par[1])+np.exp(x_plot*par3[0]+par3[1]) , 'm-', label='Fitkurve' )
plt.plot(x_plot2,np.exp(x_plot2*par3[0]+par3[1]), 'b-', label='Gerade des kurzlebigen Zerfalls' )
plt.yscale('log')
plt.xticks(
  [0, 100, 200,  400, 500, 600, 285],
  [0, 100, 200,  400, 500, 600, r'$t^*$']
)
plt.legend(loc="best")
plt.xlabel(r' $t \:/\:s$')
plt.ylabel(r' $\ln(I/Imp/15s)$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_Rh.pdf')
plt.close()

par_Rh_lang = unp.uarray(par,err1)
par_Rh_kurz = unp.uarray(par3, err3)



# tex files ##############################################################################

# tex file for a of Vanadium

with open('build/a_V.tex', 'w') as f: 
  f.write(make_SI(par_V[0]*1e03 ,r'\per\second ',exp='e-03' ,figures=2))

with open('build/b_V.tex', 'w') as f: 
  f.write(make_SI(par_V[1] ,r' ', figures=2))

with open('build/lambda_V.tex', 'w') as f: 
  f.write(make_SI(-par_V[0]*1e03 ,r'\per\second ',exp='e-03' ,figures=2))

with open('build/T_V.tex', 'w') as f: 
  f.write(make_SI(-np.log(2)/par_V[0] ,r'\second ', figures=2))

T_V = -np.log(2)/par_V[0]

with open('build/relerr_T_V.tex', 'w') as f: 
  f.write(make_SI(rel_err(T_V, 225.6) ,r'\percent ', figures=2))

# tex file for Rhodium

with open('build/a_Rh_lang.tex', 'w') as f: 
  f.write(make_SI(par_Rh_lang[0]*1e03 ,r'\per\second ', exp='e-03', figures=2))

with open('build/b_Rh_lang.tex', 'w') as f: 
  f.write(make_SI(par_Rh_lang[1] ,r' ', figures=2))

with open('build/lambda_Rh_lang.tex', 'w') as f: 
  f.write(make_SI(-par_Rh_lang[0]*1e03 ,r'\per\second ', exp='e-03', figures=2))

with open('build/T_Rh_lang.tex', 'w') as f: 
  f.write(make_SI(-np.log(2)/par_Rh_lang[0] ,r'\second ', figures=2))

T_Rh_lang =-np.log(2)/par_Rh_lang[0]

with open('build/relerr_T_Rh_lang.tex', 'w') as f: 
  f.write(make_SI(rel_err(T_Rh_lang, 258) ,r'\percent ', figures=2))


with open('build/a_Rh_kurz.tex', 'w') as f: 
  f.write(make_SI(par_Rh_kurz[0]*1e02 ,r'\per\second ', exp='e-02', figures=2))

with open('build/b_Rh_kurz.tex', 'w') as f: 
  f.write(make_SI(par_Rh_kurz[1] ,r' ', figures=2))

with open('build/lambda_Rh_kurz.tex', 'w') as f: 
  f.write(make_SI(-par_Rh_kurz[0]*1e02 ,r'\per\second ', exp='e-02', figures=2))

with open('build/T_Rh_kurz.tex', 'w') as f: 
  f.write(make_SI(-np.log(2)/par_Rh_kurz[0] ,r'\second ', figures=2))

T_Rh_kurz= -np.log(2)/par_Rh_kurz[0]

with open('build/relerr_T_Rh_kurz.tex', 'w') as f: 
  f.write(make_SI(rel_err(T_Rh_kurz, 42) ,r'\percent ', figures=2))


# Tabellen ########################################################################################

# Tabelle Für Vanadium %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1, t2 = np.array_split(t_V, 2)
N1, N2 = np.array_split(noms(N_V_err), 2)
Nerr1, Nerr2 = np.array_split(stds(N_V_err), 2)

table_header = r'''
  \begin{longtable}[H]{S[table-format=3.0] S[table-format=3.0]@{${}\pm{}$} S[table-format=2.1] S[table-format=4.0] S[table-format=2.0]@{${}\pm{}$} S[table-format=1.1] }
    \caption{Messwerte aus der Messung der Zerfälle von Vanadium. Bei den Messungen wurden
    die Zerfälle in einem Zeitintervall $\Delta t = \SI{30}{\second} $ gemessen. Die Fehler
    für die Zerfälle kommen von der Poissonverteilung und betragen $\Delta N =\sqrt{N}$. 
    }\\
    \toprule
    \multicolumn{1}{c}{ $t\:/\:s$ } & \multicolumn{2}{c}{$N\:/\:Imp/30s$ }   &
    \multicolumn{1}{c}{ $t\:/\:s$ } & \multicolumn{2}{c}{  $N\:/\:Imp/30s$ }   \\
    \cmidrule(lr{0,5em}){1-3} \cmidrule(lr{0,5em}){4-6}  

'''
table_footer = r'''    \bottomrule
  \label{tab:1}
  \end{longtable}
'''
row_template = r'     {0:1.0f} & {1:1.0f} & {2:1.1f} & {3:1.0f} & {4:1.0f} & {5:1.1f}  \\'

with open('build/V_table.tex', 'w') as g:
    g.write(table_header)
    for row in zip(t1, N1, Nerr1, t2, N2, Nerr2 ):
        g.write(row_template.format(*row).replace('nan',' ').replace('.',',') )
        g.write('\n')
    g.write(table_footer)


# Tabelle Für Rhodium %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1, t2 = np.array_split(t_Rh, 2)
N1, N2 = np.array_split(noms(N_Rh_err), 2)
Nerr1, Nerr2 = np.array_split(stds(N_Rh_err), 2)

table_header = r'''
  \begin{longtable}[H]{S[table-format=3.0] S[table-format=3.0]@{${}\pm{}$} S[table-format=2.1] S[table-format=3.0] S[table-format=2.0]@{${}\pm{}$} S[table-format=1.1] }
    \caption{Messwerte aus der Messung der Zerfälle von Rhoium. Bei den Messungen wurden
    die Zerfälle in einem Zeitintervall $\Delta t = \SI{15}{\second} $ gemessen. Die Fehler
    für die Zerfälle kommen von der Poissonverteilung und betragen $\Delta N =\sqrt{N}$. 
    }\\
    \toprule
    \multicolumn{1}{c}{ $t\:/\:s$ } & \multicolumn{2}{c}{$N\:/\:Imp/15s$ }   &
    \multicolumn{1}{c}{ $t\:/\:s$ } & \multicolumn{2}{c}{  $N\:/\:Imp/15s$ }   \\
    \cmidrule(lr{0,5em}){1-3} \cmidrule(lr{0,5em}){4-6}  

'''
table_footer = r'''    \bottomrule
  \label{tab:2}
  \end{longtable}
'''
row_template = r'     {0:1.0f} & {1:1.0f} & {2:1.1f} & {3:1.0f} & {4:1.0f} & {5:1.1f}  \\'

with open('build/Rh_table.tex', 'w') as g:
    g.write(table_header)
    for row in zip(t1, N1, Nerr1, t2, N2, Nerr2 ):
        g.write(row_template.format(*row).replace('nan',' ').replace('.',',') )
        g.write('\n')
    g.write(table_footer)
