"""Generate non cosmological population."""
import os
import numpy as np
import matplotlib.pyplot as plt
from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, unpickle, plot

MAKE = False

if MAKE:
    # Generate an FRB population
    days = 14
    population = CosmicPopulation(days*5000,
                                  z_max=0.01,
                                  lum_range=[1e40, 1e40],
                                  si_mu=0,
                                  si_sigma=0.,
                                  n_model='vol_co',
                                  pulse_model='uniform',
                                  pulse_range=[1., 1.],
                                  days=days)

    population.name = 'cosmic_z001_small'
    population.save()

    # Setup a survey
    survey = Survey('PERFECT', gain_pattern='perfect')

    # Observe the FRB population
    surv_pop = SurveyPopulation(population, survey)
    surv_pop.name = 'perfect_z001_small'
    surv_pop.save()
else:
    z4 = unpickle('perfect_z4_small')

# plot(z4, frbcat=False)
dc = np.array(z4.get('dist_co'))
z = np.array(z4.get('z'))

min_dc = min(dc)
min_z = min(z)
dc_norm = dc / min_dc
z_norm = z / min_z
plt.plot(z, dc_norm / z_norm)
plt.show()
