"""Plot options of varying FRB number densities."""
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from frbpoppy import CosmicPopulation, unpickle

MAKE = True

if MAKE:
    days = 7
    n_per_day = 5000

    # Generate population following a constant number density / comoving volume
    pop_cst = CosmicPopulation(n_per_day*days,
                               days=days,
                               z_max=6.0,
                               n_model='vol_co',
                               name='vol_co')

    # Generate population following star forming rate
    pop_sfr = CosmicPopulation(n_per_day*days,
                               days=days,
                               z_max=6.0,
                               n_model='sfr',
                               name='sfr')

    # Generate population following star forming rate
    pop_smd = CosmicPopulation(n_per_day*days,
                               days=days,
                               z_max=6.0,
                               n_model='smd',
                               name='smd')

    pop_cst.save()
    pop_sfr.save()
    pop_smd.save()

else:
    pop_cst = unpickle('vol_co')
    pop_sfr = unpickle('sfr')
    pop_sfr = unpickle('smd')


fig = plt.figure()
ax = fig.add_subplot(111)

# Get redshift of population
zs = defaultdict(list)
zs['sfr'] = pop_sfr.get('z')
zs['vol_co'] = pop_cst.get('z')
zs['smd'] = pop_smd.get('z')

for pop in zs:
    n, bins, patches = ax.hist(zs[pop],
                               50,
                               density=True,
                               label=pop,
                               histtype='step')

    plt.xlabel('Redshift')
    plt.ylabel('Probability')
    plt.legend()

# Create new legend handles but use the colors from the existing ones
handles, labels = ax.get_legend_handles_labels()
new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]

plt.legend(handles=new_handles, labels=labels)

plt.tight_layout()
plt.savefig('plots/number_density.pdf')
