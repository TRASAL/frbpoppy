"""Testing various luminosity functions."""
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
                          lum_dist_pars=[1e40, 1e50, -0.5],
                          z_max=5.0,
                          pulse=[0.1, 10],
                          si_pars=[0., 0.],
                          repeat=0.0)

    # Observe FRB population
    neg = observe(population, 'PERFECT', gain_pattern='perfect')
    neg.name = 'lum_neg'
    neg.pickle_pop()

    # Generate FRB population
    population = generate(n_per_day*days,
                          days=days,
                          lum_dist_pars=[1e45, 1e45, 0.],
                          z_max=5.0,
                          pulse=[0.1, 10],
                          si_pars=[0., 0.],
                          repeat=0.0)

    # Observe FRB population
    candle = observe(population, 'PERFECT', gain_pattern='perfect')
    candle.name = 'lum_candle'
    candle.pickle_pop()

    # Generate FRB population
    population = generate(n_per_day*days,
                          days=days,
                          lum_dist_pars=[1e40, 1e50, 0.],
                          z_max=5.0,
                          pulse=[0.1, 10],
                          si_pars=[0., 0.],
                          repeat=0.0)

    # Observe FRB population
    flat = observe(population, 'PERFECT', gain_pattern='perfect')
    flat.name = 'lum_flat'
    flat.pickle_pop()

else:
    neg = unpickle('lum_neg')
    candle = unpickle('lum_candle')
    flat = unpickle('lum_flat')


# Plot populations
plot(neg, candle, flat, mute=False, frbcat=False)
