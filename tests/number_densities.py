"""Plot options of varying FRB number densities."""
import matplotlib.pyplot as plt
import numpy as np
import os

from frbpoppy import CosmicPopulation, unpickle

MAKE = True
NUM_FRBS = True
DENSITY_FRBS = True

pop = {}
pop_types = ('cst', 'sfr', 'smd', 'shallow', 'steep')
titles = ('Constant', 'SFR', 'SMD', r'$\alpha_{in}=-0.5$',
          r'$\alpha_{in}=-2.0$')

if MAKE:
    n_frbs = int(1e6)

    # Generate population following a constant number density / comoving volume
    pop['cst'] = CosmicPopulation(n_frbs,
                                  z_max=3.,
                                  n_model='vol_co',
                                  name='cst')

    # Generate population following star forming rate
    pop['sfr'] = CosmicPopulation(n_frbs,
                                  z_max=3.,
                                  n_model='sfr',
                                  name='sfr')

    # Generate population following stellar mass density
    pop['smd'] = CosmicPopulation(n_frbs,
                                  z_max=3.,
                                  n_model='smd',
                                  name='smd')

    # Generate population following stellar mass density
    pop['shallow'] = CosmicPopulation(n_frbs,
                                      z_max=3.,
                                      n_model='vol_co',
                                      alpha=-0.5,
                                      name='shallow')

    # Generate population following stellar mass density
    pop['steep'] = CosmicPopulation(n_frbs,
                                    z_max=3.,
                                    n_model='vol_co',
                                    alpha=-2.0,
                                    name='steep')

    for k, v in pop.items():
        v.save()

else:
    for s in pop_types:
        pop[s] = unpickle(s)

# Change working directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Use A&A styling for plots
plt.style.use('./aa.mplstyle')

if NUM_FRBS:
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Get redshift of population
    zs = {}
    ns = {}
    i = 0
    for s in pop_types:
        zs[s] = pop[s].frbs.z
        ns[s], bins = np.histogram(zs[s], bins=50)

        # Remove lowest bins due to noise
        ns[s] = ns[s][1:]
        bins = bins[1:]

        bincentres = (bins[:-1] + bins[1:]) / 2

        title = titles[i]
        plt.step(bincentres, ns[s], where='mid', label=title)

        i += 1

    plt.xlabel('$z$')
    plt.ylabel(r'$\text{d}n_{\text{FRB}}/\text{d}z$')
    plt.yscale('log')
    plt.legend()
    plt.tight_layout()
    plt.savefig('plots/number_frbs.pdf')
    plt.clf()

if DENSITY_FRBS:
    plt.style.use('./aa.mplstyle')  # Use A&A styling for plots

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Get dV
    import frbpoppy.precalc as pc
    d = pc.DistanceTable().lookup(z=bincentres)
    dvols = d[3]
    dens = {}
    i = 0
    for s in pop_types:
        dens[s] = ns[s] / dvols

        title = titles[i]
        plt.step(bincentres, dens[s]/dens[s][0], where='mid', label=title)
        i += 1

    plt.xlabel('$z$')
    plt.ylabel(r'$\rho_{\text{FRB}} / \rho_{\text{FRB}}(0)$')
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig('plots/density_frbs.pdf')
