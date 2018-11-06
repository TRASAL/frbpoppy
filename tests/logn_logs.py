"""Check the log N log S slope of a population."""
from frbpoppy.population import unpickle

MAKE = True

if MAKE:
    from frbpoppy import CosmicPopulation, Survey, SurveyPopulation

    # Generate an FRB population
    days = 2
    population = CosmicPopulation(days*5000, z_max=0.5, days=days, name='example')

    # Setup a survey
    survey = Survey('PERFECT', gain_pattern='perfect')

    # Observe the FRB population
    surv_pop = SurveyPopulation(population, survey)
    surv_pop.name = 'lognlogs'
    surv_pop.save()

else:
    surv_pop = unpickle('lognlogs')

alpha, alpha_err, norm = surv_pop.calc_logn_logs()
print(alpha)
