import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
from uncertainties.unumpy import (
    nominal_values as noms,
    std_devs as stds
)

x, F = np.genfromtxt('data.txt', unpack=True)

D = F / x



table_header = r'''
  \begin{tabular}{c c c}
    \toprule
    {$\Delta x \:/\: \si{\centi\meter}$} & {$F \:/\: \si{\newton}$} & {$D \:/\: \si{\newton\per\centi\meter\tothe{-1}}$}\\
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