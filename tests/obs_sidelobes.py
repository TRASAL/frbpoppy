"""Compare detections for various numbers of sidelobes."""
from frbpoppy.do_populate import generate
from frbpoppy.do_survey import observe
from frbpoppy.do_plot import plot
from frbpoppy.population import unpickle
from frbpoppy.log import pprint

MAKE = True
MAX_SIDELOBES = 6

pops = []

if MAKE:
    days = 14
    n_per_day = 5000

    # Generate FRB population
    population = generate(n_per_day*days,
                          days=days,
                          lum_dist_pars=[1e45, 1e45, 0.],
                          z_max=2.5,
                          pulse=[0.1, 10],
                          si_pars=[0., 0.],
                          repeat=0.0)

    # Observe FRB population
    for i in range(MAX_SIDELOBES):
        pprint(f'Detecting frbs with {i} sidelobes')
        airy = observe(population, 'HTRU', gain_pattern='airy', sidelobes=i)
        airy.name = f'airy_sidelobe{i}'
        airy.pickle_pop()
        pops.append(airy)

else:
    for i in range(MAX_SIDELOBES):
        pops.append(unpickle(f'airy_sidelobe{i}'))

# Plot populations
plot(*pops, mute=False, frbcat=False)
