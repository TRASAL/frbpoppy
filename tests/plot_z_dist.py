"""Plot the distribution of redshifts following the star formation rate."""
from matplotlib import pyplot as plt

from frbpoppy.distributions import z_from_sfr

zs = []
for i in range(100000):
    zs.append(z_from_sfr(z_max=6))

plt.hist(zs)
plt.show()
