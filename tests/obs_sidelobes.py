"""Compare detections for various numbers of sidelobes."""
from frbpoppy.do_populate import generate
from frbpoppy.do_survey import observe
from frbpoppy.do_plot import plot
from frbpoppy.population import unpickle
from frbpoppy.log import pprint

OBSERVE = True
USE_SAVED = True
SIDELOBES = [0, 4, 8]

pops = []

if OBSERVE:

    if USE_SAVED:
        population = unpickle('medium')
        population.name = 'medium'
    else:
        days = 3
        n_per_day = 5000

        # Generate FRB population
        population = generate(n_per_day*days,
                              days=days,
                              lum_dist_pars=[1e45, 1e45, 0.],
                              z_max=5.0,
                              pulse=[1, 1],
                              si_pars=[0., 0.],
                              repeat=0.0)
        population.name = 'medium'
        population.pickle_pop()

    # Observe FRB population
    for i in SIDELOBES:
        pprint(f'Detecting frbs with {i} sidelobes')
        airy = observe(population,
                       'PERFECT_SMALL',
                       gain_pattern='airy',
                       sidelobes=i,
                       equal_area=SIDELOBES[-1])
        airy.name = f'airy_sidelobe{i}'
        airy.pickle_pop()
        pops.append(airy)

else:
    for i in SIDELOBES:
        pops.append(unpickle(f'airy_sidelobe{i}'))

# Plot populations
plot(*pops, mute=False, frbcat=False)
