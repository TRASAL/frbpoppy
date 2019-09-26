"""Check a powerlaw implementation generates the right numbers."""
import numpy as np
import matplotlib.pyplot as plt
import os

from frbpoppy import distributions

POWER = -1

pl = distributions.powerlaw(1e35, 1e40, POWER, 100000)
minx = min(pl)
maxx = max(pl)
bins = 10 ** np.linspace(np.log10(minx), np.log10(maxx), 50)
hist, edges = np.histogram(pl, bins=bins)

# Change working directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Use A&A styling for plots
plt.style.use('./aa.mplstyle')

fig = plt.figure()
ax = fig.add_subplot(111)
plt.step(edges[:-1], hist, where='post')
plt.xscale('log')
plt.yscale('log')
plt.savefig('plots/powerlaw.pdf')
