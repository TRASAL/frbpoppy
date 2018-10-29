"""Plot options of varying FRB number densities."""
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from frbpoppy import CosmicPopulation, unpickle

MAKE = False

if MAKE:
    days = 28
    n_per_day = 5000

    # Generate population following a constant number density / comoving volume
    pop_cst = CosmicPopulation(n_per_day*days,
                               days=days,
                               z_max=3.0,
                               n_model='constant',
                               name='constant')

    # Generate population following star forming rate
    pop_sfr = CosmicPopulation(n_per_day*days,
                               days=days,
                               z_max=3.0,
                               n_model='sfr',
                               name='sfr')

    # Generate population following star mass density
    pop_smd = CosmicPopulation(n_per_day*days,
                               days=days,
                               z_max=6.0,
                               n_model='smd',
                               name='smd')

    pop_cst.save()
    pop_sfr.save()
    pop_smd.save()

else:
    pop_cst = unpickle('constant')
    pop_sfr = unpickle('sfr')
    pop_smd = unpickle('smd')


fig = plt.figure()
ax = fig.add_subplot(111)

# Get redshift of population
zs = {}
zs['sfr'] = np.array(pop_sfr.get('z'))
zs['smd'] = np.array(pop_smd.get('z'))
zs['constant'] = np.array(pop_cst.get('z'))

n_sfr, bins = np.histogram(zs['sfr'], bins=50, density=False)
n_smd, bins = np.histogram(zs['smd'], bins=50, density=False)
n_constant, bins = np.histogram(zs['constant'], bins=50, density=False)

bincentres = [(bins[i]+bins[i+1])/2. for i in range(len(bins)-1)]

plt.step(bincentres, n_sfr/n_sfr[0], where='mid', label='SFR model')
plt.step(bincentres, n_smd/n_smd[0], where='mid', label='SMD model')
plt.step(bincentres, n_constant/n_constant, where='mid', label='Constant model')

plt.xlabel('Redshift')
plt.ylabel(r'$\rho_{FRB}(z)/\rho_{FRB}(0)$')
plt.yscale('log')
# plt.ylim(1e-1, 1e1)
plt.legend()

# Create new legend handles but use the colors from the existing ones
# handles, labels = ax.get_legend_handles_labels()
# new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]

# plt.legend(handles=new_handles, labels=labels)

plt.tight_layout()
plt.savefig('plots/number_densities.pdf')
