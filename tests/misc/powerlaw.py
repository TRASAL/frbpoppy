"""Check a powerlaw implementation generates the right numbers."""
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import gen_dists

from tests.convenience import plot_aa_style, rel_path

POWER = -1
SIZE = 1e5

pl = gen_dists.powerlaw(1e35, 1e40, POWER, int(SIZE))
minx = min(pl)
maxx = max(pl)
bins = 10 ** np.linspace(np.log10(minx), np.log10(maxx), 50)
hist, edges = np.histogram(pl, bins=bins)

plot_aa_style()

fig = plt.figure()
ax = fig.add_subplot(111)
plt.step(edges[:-1], hist, where='mid')

# Overplot the expected line
xlims = ax.get_xlim()
xs = np.array([np.log10(edges[0]), np.log10(edges[-2])])
ys = POWER*(xs-xs[0])+np.log10(hist[0])
ax.plot(10**xs, 10**ys, 'r:', label=f'Expected distribution')

plt.legend()

plt.xscale('log')
plt.yscale('log')
ax.set_ylim(1, SIZE)

plt.grid()
plt.savefig(rel_path('plots/powerlaw.pdf'))
