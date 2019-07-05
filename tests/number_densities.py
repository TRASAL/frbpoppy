"""Plot options of varying FRB number densities."""
import matplotlib.pyplot as plt
import numpy as np

from frbpoppy import CosmicPopulation, unpickle

MAKE = False
NUM_FRBS = True
DENSITY_FRBS = True

if MAKE:
    days = 1
    n_per_day = int(1e6)

    # Generate population following a constant number density / comoving volume
    pop_cst = CosmicPopulation(n_per_day*days,
                               days=days,
                               z_max=3.,
                               n_model='vol_co',
                               name='vol_co')

    # Generate population following star forming rate
    pop_sfr = CosmicPopulation(n_per_day*days,
                               days=days,
                               z_max=3.,
                               n_model='sfr',
                               name='sfr')

    # Generate population following stellar mass density
    pop_smd = CosmicPopulation(n_per_day*days,
                               days=days,
                               z_max=3.,
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

# Remove lowest bins due to noise
n_sfr = n_sfr[1:]
n_constant = n_constant[1:]
n_smd = n_smd[1:]
bins = bins[1:]

bincentres = (bins[:-1] + bins[1:]) / 2

if NUM_FRBS:
    plt.step(bincentres, n_sfr, where='mid', label='SFR')
    plt.step(bincentres, n_constant, where='mid', label='Constant')
    plt.step(bincentres, n_smd, where='mid', label='SMD')

    plt.xlabel('$z$')
    plt.ylabel(r'$n_{\text{FRB}}$')
    plt.yscale('log')
    plt.legend()
    plt.tight_layout()
    plt.savefig('plots/number_frbs.pdf')
    plt.clf()

if DENSITY_FRBS:
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Get dV
    import frbpoppy.precalc as pc
    d = pc.DistanceTable().lookup(z=bincentres)
    dvols = d[3]
    den_sfr = n_sfr / dvols
    den_con = n_constant / dvols
    den_smd = n_smd / dvols

    plt.step(bincentres, den_sfr/den_sfr[0], where='mid', label='SFR')
    plt.step(bincentres, den_con/den_con[0], where='mid', label='Constant')
    plt.step(bincentres, den_smd/den_smd[0], where='mid', label='SMD')

    plt.xlabel('$z$')
    plt.ylabel(r'$\rho_{\text{FRB}} / \rho_{\text{FRB}}(0)$')
    plt.yscale('log')
    plt.legend()
    plt.tight_layout()
    plt.savefig('plots/density_frbs.pdf')
