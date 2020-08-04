"""Check a powerlaw implementation generates the right numbers."""
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import gen_dists

from tests.convenience import plot_aa_style, rel_path

POWER = -1

pl = gen_dists.powerlaw(1e35, 1e40, POWER, 100000)
minx = min(pl)
maxx = max(pl)
bins = 10 ** np.linspace(np.log10(minx), np.log10(maxx), 50)
hist, edges = np.histogram(pl, bins=bins)

plot_aa_style()

fig = plt.figure()
ax = fig.add_subplot(111)
plt.step(edges[:-1], hist, where='post')
plt.xscale('log')
plt.yscale('log')
plt.savefig(rel_path('plots/powerlaw.pdf'))
