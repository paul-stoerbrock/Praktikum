import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
from uncertainties.unumpy import (
    nominal_values as noms,
    std_devs as stds
)


def linreg(x, y):

    m= (np.mean(x*y)-np.mean(x)* np.mean(y))/(np.mean(x^2)-np.mean(x)^2)
    b= (np.mean(y)* np.mean(x^2)-np.mean(xy)* np.mean(x))/(np.mean(x^2)-np.mean(x)^2)

    return m,b

x, F = np.genfromtxt('data.txt', unpack=True)

D = F / x

params, cov_mat = np.polyfit(x, F, deg=1, cov=True)
plt.plot(x, F, 'k.', label="Messwerte")

x_plot = np.linspace(0, 10)
plt.plot(
    x_plot,
    params[0] * x_plot + params[1],
    label ='Lineare Regression'
)
plt.legend(loc="best")

plt.savefig('build/plot.pdf')

D_bar=np.mean(D)

meanD=r' $\SI{{D_bar}}{\centi\meter}$'

with open('build/meanD.tex', 'w') as g:
    g.write(meanD)







table_header = r'''
  \begin{tabular}{c c c}
    \toprule
    {$\Delta x \:/\: \si{\centi\meter}$} & {$F \:/\: \si{\newton}$} & {$D \:/\: \si{\newton\per\centi\meter}$}\\
    \midrule
'''

table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'     {0:1.1f} & {1:1.2f} & {2:1.2f}  \\'


with open('build/table.tex', 'w') as f:
    f.write(table_header)
    for row in zip(x, F, D):
        f.write(row_template.format(*row))
        f.write('\n')
    f.write(table_footer)