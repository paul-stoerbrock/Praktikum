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
import math
import sympy

def make_SI(num, unit, exp='', figures=None):
    ''' Format an uncertainties ufloat as a \SI quantity '''
    if np.any(stds([num])):
        if figures is None:
            figures = ''
        x = '{0:.{1:}uf}'.format(num, figures).replace('/', '')
    else:
        x = '{0:.{1:}f}'.format(num, figures)

    return r'\SI{{{}{}}}{{{}}}'.format(x, exp, unit)


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
        
    return latex(sympy.sqrt(s), symbol_names=latex_names) # automatic Gauß error function

def rel_err(mess, lit):
    return abs((unp.nominal_values(mess)-unp.nominal_values(lit))/unp.nominal_values(lit))*100

def m( U, I):
    m_T = np.array([])
    x=0
    while x < len(U)-1:
        m_T = np.append(m_T, (I[x+1]-I[x])/(U[x+1]-U[x]))
        x=x+1    
    return m_T      # Funktion zur automatischen Berechnung der Steigung für die differentielle Energieverteilung der Elektronen

def p(T):
    return 5.5*10**7 *np.exp(-6876/T)

def w(T):
    return 0.0029/p(T)

# Messwerte ######################################################################

UA_T1, IA_UA_T1 = np.genfromtxt('dataUA_T1.txt' ,unpack = True) # UA_T1 bei 25,3°C gemessen bei UB =11V, UA_T2 bei 145,5°C
UA_T2, IA_UA_T2 = np.genfromtxt('dataUA_T2.txt' ,unpack = True) # UA_T2 bei 145,5°C
UB_T1 , IA_UB_T1, UB_T2, IA_UB_T2 = np.genfromtxt('dataUB.txt' ,unpack = True) # UB_T1 bei 164°C und UA=1V , UB_T2 bei 175°C
 
T1_UA = const.convert_temperature( 25.3, 'C' ,'K')
T2_UA = const.convert_temperature( 145.5, 'C' ,'K')
U_B = 11

T1_UB = const.convert_temperature( 164, 'C' ,'K')
T2_UB = const.convert_temperature( 175, 'C' ,'K')
U_A = 1
print(T1_UB)
print(T2_UB)

m_T1 = m(UA_T1, IA_UA_T1)
m_T2 = m(UA_T2[:5], IA_UA_T2[:5])

# Berechnung des Wirkungsquerschnitt

T = np.array([T1_UA, T2_UA, T1_UB, T2_UB])
p_sat = p(T)
w_bar = w(T)


# Plots ####################################################################


# Plot für die Steigung der Gegenspannung %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



plt.plot(UA_T1[:10], m_T1 ,'kx' , label='Messwerte bei 298,4 K')
plt.legend(loc="best")
plt.xlabel(r' $U_A \:/\:\mathrm{V}$')
plt.ylabel(r'$m \:/\: \mathrm{A/V}$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_UA.pdf')
plt.close()


plt.plot(UA_T2[:4], m_T2 ,'kx' , label='Messwerte bei 418,6 K')
plt.legend(loc="best")
plt.xlabel(r' $U_A \:/\:\mathrm{V}$')
plt.ylabel(r'$m \:/\: \mathrm{A/V}$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_UA_2.pdf')
plt.close()






# Plot für Die Beschleunigungsspannung %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plt.plot(UB_T1, IA_UB_T1 ,'kx' , label='Messwerte bei 437,1 K')
plt.plot(UB_T2, IA_UB_T2 ,'bx' , label='Messwerte bei 448,1 K')
plt.axvline(52, c='g',label='Maxima bei 448,1 K')
plt.axvline(30, c='g')
plt.axvline(36, c='g')
plt.axvline(41, c='g',)
plt.axvline(46, c='g')
plt.legend(loc="best")
plt.xlabel(r' $U_B \:/\:\mathrm{V}$')
plt.ylabel(r'$I_A \:/\: \mathrm{nA}$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_UB.pdf')
plt.close()


U1 = np.array([6,5,5,6])
U1_err = ufloat(np.mean(U1), sem(U1))




# Tex Files ####################################################################


# tex file of U1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

with open('build/U1.tex', 'w') as f:
  f.write(make_SI(U1_err ,r'\electronvolt', figures=2))

# tex file of relative error U1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

with open('build/relerr_U1.tex', 'w') as f:
  f.write(make_SI(rel_err(U1_err, 4.9) ,r'\percent', figures=2))

# tex file of lambda %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

with open('build/lambda.tex', 'w') as f:
  f.write(make_SI(const.h *const.c/(U1_err*const.e) *1e9 ,r'\nano\meter', figures=2))

# Tabellen ######################################################################


# Tabelle der Temperaturen sowie Druck und Wirkungsquerschnitt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table_header = r'''
  \begin{longtable}[H]{S[table-format=3.2] S[table-format=5.1] S[table-format=3.1] S[table-format=3.1] }
    \caption{
        In der Tabelle sind die Temperaturen mitsamt ihres Sättigungsdampfdruck und des
        Wirkungsquerschnitts dargestellt. Der Sättigungsdampfdruck wird mithilfe der Formel
        \eqref{eq:1} und der Wirkungsquerschnitt mit Formel \eqref{eq:5} berechnet. In der vierten
        Spalte ist dann das Verhältnis zwischen dem Abstand des Glühdrahts zur Beschleunigerelektrode a,
        welche hier 1 cm beträgt, gegen den Wirkungsquerschnitt aufgetragen.
    }\\ 
    \toprule
    \multicolumn{1}{c}{ $T\:/\:K$ } & \multicolumn{1}{c}{$p_{\text{sät}} \:/\:\mu bar$ }   &
    \multicolumn{1}{c}{ $\bar w \:/\:\mu m$ } &  \multicolumn{1}{c}{ $\frac{a}{\bar w} $ }   \\
    \cmidrule(lr{0,5em}){1-4} 

'''
table_footer = r'''    \bottomrule
  \label{tab:1}
  \end{longtable}
'''
row_template = r'     {0:1.1f} & {1:1.2f} & {2:1.1f} & {3:3.1f}   \\'

with open('build/T_table.tex', 'w') as g:
    g.write(table_header)
    for row in zip(T, p_sat*1e03, w_bar*1e04, 1/w_bar):
        g.write(row_template.format(*row).replace('nan',' ').replace('.',',') )
        g.write('\n')
    g.write(table_footer)


# Tabelle der Gegenspannung %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table_header = r'''
  \begin{longtable}[H]{S[table-format=2.0] S[table-format=3.0] S[table-format=2.0] S[table-format=3.0] }
    \caption{
        Messwerte für die differentielle Energieverteilung.
        Gemessen wird die Gegenspannung $U_A$ gegen den Auffängerstrom
        $I_A$ bei den Temperaturen \SI{298.4}{\kelvin} und \SI{418.6}{\kelvin}.
    }\\ 
    \toprule
    \multicolumn{2}{c}{T=\SI{298.4}{\kelvin}} & \multicolumn{2}{c}{T=\SI{418.6}{\kelvin}} \\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4}  \\
    \multicolumn{1}{c}{ $U_A\:/\:\si{\volt} $ } & \multicolumn{1}{c}{$I\:/\:\si{\nano\ampere} $ }   &
    \multicolumn{1}{c}{ $U_A\:/\:\si{\volt}$ } & \multicolumn{1}{c}{$I\:/\:\si{\nano\ampere}$ }    \\
    \cmidrule(lr{0,5em}){1-4}  

'''
table_footer = r'''    \bottomrule
  \label{tab:2}
  \end{longtable}
'''
row_template = r'     {0:2.2f} & {1:3.2f} & {2:2.2f} & {3:3.2f}   \\'

with open('build/UA_table.tex', 'w') as g:
    g.write(table_header)
    for row in zip( UA_T1, IA_UA_T1*1e9, UA_T2, IA_UA_T2*1e9):
        g.write(row_template.format(*row).replace('nan',' ').replace('.',',') )
        g.write('\n')
    g.write(table_footer)




# Tabelle für Frank-Hertz-Kurve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IA_UB_T1_1, IA_UB_T1_2 = np.array_split(IA_UB_T1*1e12, 2)
IA_UB_T2_1, IA_UB_T2_2 = np.array_split(IA_UB_T2*1e12, 2)
UB_T1_1, UB_T1_2 = np.array_split(UB_T1, 2)
UB_T2_1, UB_T2_2 = np.array_split(UB_T2, 2)

table_header = r'''
  \begin{longtable}[H]{S[table-format=2.0] S[table-format=3.0] S[table-format=2.0] S[table-format=3.0] S[table-format=2.0] S[table-format=3.0] S[table-format=2.0] S[table-format=3.0]  }
    \caption{
        In der Tabelle werden die Messwerte der Frank-Hertz-Kurve für die beiden
        Temperaturen \SI{437.5}{\kelvin} und \SI{448.1}{\kelvin} dargestellt.
        Gemessen wird die Beschleunigungspannung gegen den Auffängerstrom.
    }\\ 
    \toprule
    \multicolumn{2}{c}{T=\SI{437.15}{\kelvin}} & \multicolumn{2}{c}{T=\SI{448.15}{\kelvin}} & \multicolumn{2}{c}{T=\SI{437.15}{\kelvin}} & \multicolumn{2}{c}{T=\SI{448.15}{\kelvin}} \\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4} \cmidrule(lr{0,5em}){5-6} \cmidrule(lr{0,5em}){7-8} \\
    \multicolumn{1}{c}{ $U_B\:/\:\si{\volt} $ } & \multicolumn{1}{c}{$I\:/\:\si{\pico\ampere}$ }   &
    \multicolumn{1}{c}{ $U_B\:/\:\si{\volt} $ } & \multicolumn{1}{c}{$I\:/\:\si{\pico\ampere}$ }   &
    \multicolumn{1}{c}{ $U_B\:/\:\si{\volt} $ } & \multicolumn{1}{c}{$I\:/\:\si{\pico\ampere}$ }   &
    \multicolumn{1}{c}{ $U_B\:/\:\si{\volt} $ } & \multicolumn{1}{c}{$I\:/\:\si{\pico\ampere}$ }  \\
    \cmidrule(lr{0,5em}){1-4} \cmidrule(lr{0,5em}){5-8}  

'''
table_footer = r'''    \bottomrule
  \label{tab:3}
  \end{longtable}
'''
row_template = r'     {0:2.0f} & {1:3.0f} & {2:2.0f} & {3:3.0f} & {4:2.0f} & {5:3.0f} & {6:2.0f} & {7:3.0f}  \\'

with open('build/UB_table.tex', 'w') as g:
    g.write(table_header)
    for row in zip(UB_T1_1, IA_UB_T1_1, UB_T2_1, IA_UB_T2_1, UB_T1_2, IA_UB_T1_2, UB_T2_2, IA_UB_T2_2  ):
        g.write(row_template.format(*row).replace('nan',' ').replace('.',',') )
        g.write('\n')
    g.write(table_footer)


