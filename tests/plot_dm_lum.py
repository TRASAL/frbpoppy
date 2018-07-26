"""Plot opitions of varying FRB number densities."""
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

from frbpoppy.do_survey import observe
from frbpoppy.do_populate import generate
from frbpoppy.population import unpickle

MAKE = False

if MAKE:
    days = 28
    n_per_day = 5000

    # Generate FRB population
    population = generate(n_per_day*days,
                          days=days,
                          lum_range=[1e40, 1e50],
                          lum_index=-0.25,
                          z_max=2.5,
                          pulse_model='uniform',
                          pulse_range=[5., 5.],
                          si_mu=0.,
                          si_sigma=0.,
                          repeat=0.0)

    # Observe FRB population
    neg = observe(population, 'HTRU', gain_pattern='perfect')
    neg.name = 'lum_neg'
    neg.pickle_pop()

    # Generate FRB population
    population = generate(n_per_day*days,
                          days=days,
                          lum_range=[1e45, 1e45],
                          lum_index=0.,
                          z_max=2.5,
                          pulse_model='uniform',
                          pulse_range=[5., 5.],
                          si_mu=0.,
                          si_sigma=0.,
                          repeat=0.0)

    # Observe FRB population
    candle = observe(population, 'HTRU', gain_pattern='perfect')
    candle.name = 'lum_candle'
    candle.pickle_pop()

    # Generate FRB population
    population = generate(n_per_day*days,
                          days=days,
                          lum_range=[1e40, 1e50],
                          lum_index=0.,
                          z_max=2.5,
                          pulse_model='uniform',
                          pulse_range=[5., 5.],
                          si_mu=0.,
                          si_sigma=0.,
                          repeat=0.0)
    # Observe FRB population
    flat = observe(population, 'HTRU', gain_pattern='perfect')
    flat.name = 'lum_flat'
    flat.pickle_pop()

else:
    neg = unpickle('lum_neg')
    candle = unpickle('lum_candle')
    flat = unpickle('lum_flat')

neg.name = 'Negative'
candle.name = 'Candle'
flat.name = 'Flat'

# Get dispersion measure of population
dms = defaultdict(list)

for pop in (neg, candle, flat):
    for src in pop.sources:
        dms[pop.name].append(src.dm)

for pop in dms:
    dmp = dms[pop]

    # Bin up
    hist, edges = np.histogram(dmp, bins=np.arange(0, 2500, 50))
    cum_hist = [sum(hist[:i]) for i in range(len(hist))]

    # Normalise to one
    m = max(cum_hist)
    cum_hist = [c/m for c in cum_hist]

    plt.plot(edges[:-1], cum_hist, label=pop)

    plt.xlabel(r'Dispersion Measure ($\text{pc}\ \text{cm}^{-3}$)')
    plt.ylabel(r'N(<DM)')
    plt.legend()

plt.yscale('log')
plt.tight_layout()
plt.savefig('plots/dm_lum.pdf')
