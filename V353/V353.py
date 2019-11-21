
import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
from uncertainties.unumpy import (
    nominal_values as noms,
    std_devs as stds
)


x = np.linspace(0, 1)
plt.plot(x, x**2, 'b-')
plt.savefig('build/plot.pdf')