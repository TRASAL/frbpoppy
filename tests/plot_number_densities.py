"""Plot options of varying FRB number densities."""
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from frbpoppy import CosmicPopulation, unpickle

MAKE = False

if MAKE:
    days = 50
    n_per_day = 5000

    # Generate population following a constant number density / comoving volume
    pop_cst = CosmicPopulation(n_per_day*days,
                               days=days,
                               z_max=3.0,
                               n_model='vol_co',
                               name='vol_co')

    # Generate population following star forming rate
    pop_sfr = CosmicPopulation(n_per_day*days,
                               days=days,
                               z_max=3.0,
                               n_model='sfr',
                               name='sfr')

    # Generate population following stellar mass density
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
    pop_smd = unpickle('smd')


fig = plt.figure()
ax = fig.add_subplot(111)

# Get redshift of population
zs = {}
zs['sfr'] = pop_sfr.frbs.z
zs['smd'] = pop_smd.frbs.z
zs['vol_co'] = pop_cst.frbs.z

n_sfr, bins = np.histogram(zs['sfr'], bins=100, density=False)
n_smd, bins = np.histogram(zs['smd'], bins=100, density=False)
n_constant, bins = np.histogram(zs['vol_co'], bins=100, density=False)

bincentres = (bins[:-1] + bins[1:]) / 2

plt.step(bincentres, n_sfr, where='mid', label='SFR')
plt.step(bincentres, n_constant, where='mid', label='Constant')
plt.step(bincentres, n_smd, where='mid', label='SMD')

plt.xlabel('Redshift')
plt.ylabel(r'$\rho_{FRB}(z)$')
plt.yscale('log')
plt.legend()
plt.tight_layout()
plt.savefig('plots/number_densities.pdf')
