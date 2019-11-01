import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
from uncertainties.unumpy import (
    nominal_values as noms,
    std_devs as stds,

t, U, U_err = np.genfromtxt('data.txt', unpack=True)
t *= 1e-3
U = 1e3 * unp.uarray(U, U_err)

table_header = r'''
  \begin{tabular}{c c c}
    \toprule
    {$/Delta x \:/\: [cm]$} & \multicolumn{2}{c}{$F \:/\: [N]$} & \multicolumn{3}{c}{$D$}
    \midrule
'''

table_footer = r'''    \bottomrule
  \end{tabular}
'''
row_template = r'    {0:1.3f} & {1.n:1.2f} & {1.s:1.2f} & {2:1.3f} & {3.n:1.2f} & {3.s:1.2f} \\'


with open('build/loesung-table.tex', 'w') as f:
    f.write(table_header)
    for row in zip(t1, U1, t2, U2):
        f.write(row_template.format(*row))
        f.write('\n')
    f.write(table_footer)