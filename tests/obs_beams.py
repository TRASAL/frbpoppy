"""Compare detections for various beam patterns."""
from frbpoppy.do_populate import generate
from frbpoppy.do_survey import observe
from frbpoppy.do_plot import plot
from frbpoppy.population import unpickle

MAKE = False

if MAKE:
    days = 7
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
    gaussian = observe(population, 'HTRU', gain_pattern='gaussian')
    gaussian.name = 'gaussian'
    gaussian.pickle_pop()

    # Observe FRB population
    airy = observe(population, 'HTRU', gain_pattern='airy')
    airy.name = 'airy'
    airy.pickle_pop()

    # Observe FRB population
    parkes = observe(population, 'HTRU', gain_pattern='parkes')
    parkes.name = 'parkes'
    parkes.pickle_pop()

else:
    gaussian = unpickle('gaussian')
    airy = unpickle('airy')
    parkes = unpickle('parkes')

# Plot populations
plot(gaussian, airy, parkes, mute=False, frbcat=False)
