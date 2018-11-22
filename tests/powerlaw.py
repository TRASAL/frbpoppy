import numpy as np
import matplotlib.pyplot as plt
import frbpoppy.distributions as dis
from frbpoppy.number_density import NumberDensity

n = NumberDensity(model='vol_co')
s = [n.draw() for i in range(int(1e4))]
# The error happens somewhere in the vectorizing of the function
# N(S)
number, bins = np.histogram(np.log10(s), bins=500)
# N(>S) from N(S)
n_gt_s = np.cumsum(number[::-1])[::-1]
# logS
x = bins[:-1]
y = np.log10(number)


plt.plot(x, y)
plt.show()
