"""Check ASKAP predictions."""
from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, plot, unpickle

MAKE = False
d = {'askap': None, 'htru': None, 'arecibo': None}

if MAKE:
    # Generate an FRB population
    days = 7
    population = CosmicPopulation(days*5000,
                                  z_max=1.0,
                                  lum_range=[1e45, 1e45],
                                  si_mu=0,
                                  si_sigma=0.,
                                  n_model='constant',
                                  pulse_model='uniform',
                                  pulse_range=[1., 1.],
                                  days=days,
                                  dm_host_model='lognormal',
                                  dm_host_mu=800,
                                  dm_host_sigma=2.7)

    population.name = 'cosmic_askap_prediction'
    population.save()

    # Setup surveys
    askap = Survey('ASKAP-FLY', gain_pattern='airy', sidelobes=0.5)
    htru = Survey('HTRU', gain_pattern='parkes')
    arecibo = Survey('ARECIBO-SPF', gain_pattern='airy')

    d = {'askap': askap, 'htru': htru, 'arecibo': arecibo}

    # Observe the FRB populations
    for name, survey in d.items():
        surv_pop = SurveyPopulation(population, survey)
        surv_pop.name = f'{name}_askap_prediction'
        print(surv_pop.rates())
        surv_pop.save()
        d[name] = surv_pop
else:
    for name in d:
        d[name] = unpickle(f'{name}_askap_prediction')

# Plot populations
plot(*d.values(), frbcat=False, mute=False)
