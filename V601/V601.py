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
m_T2 = m(UA_T2, IA_UA_T2)

# Berechnung des Wirkungsquerschnitt

T = np.array([T1_UA, T2_UA, T1_UB, T2_UB])
p_sat = p(T)
w_bar = w(T)
print(T)
print(p_sat)
print(w_bar)

# Plots ####################################################################


# Plot für die Steigung der Gegenspannung %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#par1, cov1 = np.polyfit( UA_T1, IA_UA_T1, deg=1, cov=True )
#
#err1 = np.sqrt(np.diag(cov1))
#
#par2, cov2 = np.polyfit( UA_T2, IA_UA_T2, deg=1, cov=True )
#
#err2 = np.sqrt(np.diag(cov2))


plt.plot(UA_T1[:10], m_T1 ,'kx' , label='Messwerte bei 25,3°C')
#plt.plot(UA_T2[:5], m_T2 ,'bx' , label='Messwerte bei 145,5°C')
#x_plot = np.linspace(0, 5, 1000)
#plt.plot(x_plot, par1[0]*x_plot+par1[1], 'r-', label='Fitgerade bei 25,3°C')
#plt.plot(x_plot, par2[0]*x_plot+par2[1], 'b-', label='Fitgerade bei 145,5°C')
plt.legend(loc="best")
plt.xlabel(r' $U_A \:/\:V$')
plt.ylabel(r'$m \:/\: A/V$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_UA.pdf')
plt.close()


#par_T1 = unp.uarray(par1, err1)
#
#par_T2 = unp.uarray(par2, err2)



# Plot für Die Beschleunigungsspannung %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plt.plot(UB_T1, IA_UB_T1 ,'kx' , label='Messwerte bei 164°C')
plt.plot(UB_T2, IA_UB_T2 ,'bx' , label='Messwerte bei 175°C')
plt.axvline(25, c='g',label='Maxima')
plt.axvline(31, c='g')
plt.axvline(36, c='g')
plt.axvline(41, c='g',)
plt.axvline(46, c='g')
plt.legend(loc="best")
plt.xlabel(r' $U_B \:/\:V$')
plt.ylabel(r'$I_A \:/\: nA$')
plt.grid()
plt.tight_layout
plt.savefig('build/plot_UB.pdf')
plt.close()




# Tabellen ######################################################################


# Tabelle der Temperaturen sowie Druck und Wirkungsquerschnitt

table_header = r'''
  \begin{longtable}[H]{S[table-format=3.2] S[table-format=5.1] S[table-format=3.1] }
    \caption{
    }\\ 
    \toprule
    \multicolumn{1}{c}{ $T\:/\:K$ } & \multicolumn{1}{c}{$p_{\text{sät}} \:/\:\mu bar$ }   &
    \multicolumn{1}{c}{ $\bar w \:/\:\mu m$ }     \\
    \cmidrule(lr{0,5em}){1-3} 

'''
table_footer = r'''    \bottomrule
  \label{tab:1}
  \end{longtable}
'''
row_template = r'     {0:1.1f} & {1:1.2f} & {2:1.1f}   \\'

with open('build/T_table.tex', 'w') as g:
    g.write(table_header)
    for row in zip(T, p_sat*1e03, w_bar*1e04):
        g.write(row_template.format(*row).replace('nan',' ').replace('.',',') )
        g.write('\n')
    g.write(table_footer)



# Tabelle für Frank-hertz-Kurve

IA_UB_T1_1, IA_UB_T1_2 = np.array_split(IA_UB_T1*1e12, 2)
IA_UB_T2_1, IA_UB_T2_2 = np.array_split(IA_UB_T2*1e12, 2)
UB_T1_1, UB_T1_2 = np.array_split(UB_T1, 2)
UB_T2_1, UB_T2_2 = np.array_split(UB_T2, 2)

table_header = r'''
  \begin{longtable}[H]{S[table-format=2.0] S[table-format=3.0] S[table-format=2.0] S[table-format=3.0] S[table-format=2.0] S[table-format=3.0] S[table-format=2.0] S[table-format=3.0]  }
    \caption{
    }\\ 
    \toprule
    \multicolumn{2}{c}{T=437,15 K} & \multicolumn{2}{c}{T=448,15 K} & \multicolumn{2}{c}{T=437,15 K} & \multicolumn{2}{c}{T=448,15 K} \\
    \cmidrule(lr{0,5em}){1-2} \cmidrule(lr{0,5em}){3-4} \cmidrule(lr{0,5em}){5-6} \cmidrule(lr{0,5em}){7-8} \\
    \multicolumn{1}{c}{ $U_B\:/\:V$ } & \multicolumn{1}{c}{$I\:/\:pA$ }   &
    \multicolumn{1}{c}{ $U_B\:/\:V$ } & \multicolumn{1}{c}{$I\:/\:pA$ }   &
    \multicolumn{1}{c}{ $U_B\:/\:V$ } & \multicolumn{1}{c}{$I\:/\:pA$ }   &
    \multicolumn{1}{c}{ $U_B\:/\:V$ } & \multicolumn{1}{c}{$I\:/\:pA$ }  \\
    \cmidrule(lr{0,5em}){1-4} \cmidrule(lr{0,5em}){5-8}  

'''
table_footer = r'''    \bottomrule
  \label{tab:2}
  \end{longtable}
'''
row_template = r'     {0:2.0f} & {1:3.0f} & {2:2.0f} & {3:3.0f} & {4:2.0f} & {5:3.0f} & {6:2.0f} & {7:3.0f}  \\'

with open('build/UB_table.tex', 'w') as g:
    g.write(table_header)
    for row in zip(UB_T1_1, IA_UB_T1_1, UB_T2_1, IA_UB_T2_1, UB_T1_2, IA_UB_T1_2, UB_T2_2, IA_UB_T2_2  ):
        g.write(row_template.format(*row).replace('nan',' ').replace('.',',') )
        g.write('\n')
    g.write(table_footer)


