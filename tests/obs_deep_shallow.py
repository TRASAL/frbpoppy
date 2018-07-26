"""Compare detections a deep and shallow survey."""
from frbpoppy.do_populate import generate
from frbpoppy.do_survey import observe
from frbpoppy.do_plot import plot
from frbpoppy.population import unpickle

MAKE = True

if MAKE:
    days = 7
    n_per_day = 5000

    # Generate FRB population
    population = generate(n_per_day*days,
                          days=days,
                          lum_range=[1e40, 1e45],
                          lum_index=0,
                          z_max=2.5,
                          si_mean=0.,
                          si_sigma=0.)

    # Observe FRB population
    deep = observe(population, 'TEST-DEEP', gain_pattern='perfect')
    deep.name = 'deep'
    deep.pickle_pop()

    # Observe FRB population
    shallow = observe(population, 'TEST-SHALLOW', gain_pattern='perfect')
    shallow.name = 'shallow'
    shallow.pickle_pop()


else:
    deep = unpickle('deep')
    shallow = unpickle('shallow')

# Plot populations
plot(deep, shallow, mute=False, frbcat=False)
