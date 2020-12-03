"""Visualise run."""
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

from simulations import SimulationOverview
from goodness_of_fit import GoodnessOfFit

RUN = 5
run_pars = {1: ['alpha', 'si', 'li'],
            2: ['lum_min', 'lum_max', 'li'],
            3: ['w_mean', 'w_std'],
            4: ['dm_igm_slope', 'dm_host']}

so = SimulationOverview()
df = so.df[so.df.run == RUN].copy()
par_set = df.par_set.iloc[0]
cols = run_pars[par_set]

# Group together survey results
gofs = {c: [] for c in cols}
gofs['gof'] = []
gf = GoodnessOfFit()
for vals, group in df.groupby(cols):
    gofs['gof'].append(gf.weighted_median(group))
    for i, c in enumerate(cols):
        gofs[c].append(vals[i])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Log business
labels = []
for col in cols:
    diff = np.diff(np.sort(list(set(gofs[col]))))
    if not np.isclose(diff[0], diff[-1]):
        gofs[col] = np.log10(gofs[col])
        label = f'log {col}'
    else:
        label = col
    labels.append(label)

x = gofs[cols[0]]
y = gofs[cols[1]]
colours = np.array(gofs['gof'])
if len(cols) == 2:
    z = 0
else:
    z = gofs[cols[2]]
    ax.set_zlabel(labels[2])

p = ax.scatter(x, y, z, c=colours, s=20*10**colours)
cbar = fig.colorbar(p, orientation='horizontal',
                    cax=fig.add_axes([0.05, 0.9, 0.2, 0.03]))
cbar.ax.set_ylabel('GoF')  #, rotation=-90, va="bottom")

xlabel = labels[0]
ylabel = labels[1]

# ugly, but quick
if xlabel == 'alpha':
    xlabel = r'$\alpha$'
if xlabel == 'log lum_min':
    xlabel = r'lum$_{\rm min}$'
if xlabel == 'log lum_max':
    xlabel = r'lum$_{\rm max}$'
if ylabel == 'log lum_min':
    ylabel = r'lum$_{\rm min}$'
if ylabel == 'log lum_max':
    ylabel = r'lum$_{\rm max}$'
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)

plt.tight_layout()
plt.show()
