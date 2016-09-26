import math
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from numpy import trapz

# question 4
h = (1.0, 0.5, 0.1, 0.05, 0.01)
N = (8.0, 16.0, 80.0, 160.0, 800.0)
rel_err = (8.44660179e-3, 2.30286448e-3, 9.84273963e-5, 2.48043656e-5, 9.98488163e-7)

lin_fit_h = np.polyfit(np.log10(h), np.log10(rel_err), 1)
lin_fit_N = np.polyfit(np.log10(N), np.log10(rel_err), 1)

fig, ax1 = plt.subplots()
ax1.plot(np.log10(h), np.log10(rel_err), 'k*', label='Mesh spacing, m=%.2f' %lin_fit_h[0])
ax1.plot(np.log10(N), np.log10(rel_err), 'r*', label='Cell count, m=%.2f' %lin_fit_N[0])
ax1.plot(np.log10(h), lin_fit_h[0] * np.log10(h) + lin_fit_h[1], 'k--')
ax1.plot(np.log10(N), lin_fit_N[0] * np.log10(N) + lin_fit_N[1], 'r--')

legend = ax1.legend()
plt.xlabel('$log_{10}(h)$ or $log_{10}(N)$')
plt.ylabel('$log_{10}$(Relative error)')
plt.savefig("Question4.pdf", format='pdf', bbox_inches='tight')

plt.close()
