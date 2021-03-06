import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
from uncertainties.unumpy import (
    nominal_values as noms,
    std_devs as stds
)


def linreg(x, y):

    m= (np.mean(x*y)-np.mean(x)* np.mean(y))/(np.mean(x*x)-np.mean(x)*np.mean(x))
    b= (np.mean(y)* np.mean(x*x)-np.mean(x*y)* np.mean(x))/(np.mean(x*x)-np.mean(x)*np.mean(x))

    return m,b

x, F = np.genfromtxt('data.txt', unpack=True)

x *= 1e-02

D = F / x

params, cov_mat = np.polyfit(x, F, deg=1, cov=True)
plt.plot(x, F, 'k.', label="Messwerte")

x_plot = np.linspace(0, 60)
plt.plot(
    x_plot,
    params[0] * x_plot + params[1],
    'r',
    label ='Lineare Regression'
)
plt.legend(loc="best")
plt.title('Linear Regression')
plt.xlabel('Auslenkung in cm')
plt.ylabel('Newton')
plt.savefig('build/plot.pdf')

D_bar=np.mean(D)

meanD=r' }{\newton\per\meter}$'

with open('build/meanD.tex', 'w') as g:
    g.write('$\SI{')
    g.write(f'{D_bar:.2f}')
    g.write(meanD)

D_linreg, b=linreg(x, F)


linregD=r' }{\newton\per\meter}$'

with open('build/linregD.tex','w') as h:
    h.write('$\SI{')
    h.write(f'{D_linreg:.2f}')
    h.write(linregD)


table_header = r'''
  \begin{tabular}{c c c}
    \toprule
    {$\Delta x \:/\: \si{\meter}$} & {$F \:/\: \si{\newton}$} & {$D \:/\: \si{\newton\meter\tothe{-1}}$}\\
    \midrule
'''

table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.2f} & {1:1.2f} & {2:1.2f}  \\'


with open('build/table.tex', 'w') as f:
    f.write(table_header)
    for row in zip(x, F, D):
        f.write(row_template.format(*row))
        f.write('\n')
    f.write(table_footer)