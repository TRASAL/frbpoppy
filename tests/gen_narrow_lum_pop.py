"""Generate non cosmological population."""
import os
from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, unpickle

MAKE = False

if MAKE:
    # Generate an FRB population
    days = 30
    population = CosmicPopulation(days*5000,
                                  z_max=0.01,
                                  lum_range=[1e39, 1e40],
                                  si_mu=0,
                                  si_sigma=0.,
                                  n_model='constant',
                                  pulse_range=[1., 1.],
                                  days=days)

    population.name = 'cosmic_z001_lum_narrow'
    population.save()

    # Setup a survey
    survey = Survey('PERFECT', gain_pattern='perfect')

    # Observe the FRB population
    surv_pop = SurveyPopulation(population, survey)
    surv_pop.name = 'perfect_z001_lum_narrow'
    surv_pop.save()
else:
    surv_pop = unpickle('perfect_z001_lum_narrow')

fluences = surv_pop.get('fluence')

with open(os.path.expanduser('~/Downloads/fluence_z001.csv'), 'w') as f:
    for item in fluences:
        f.write("%s\n" % item)
