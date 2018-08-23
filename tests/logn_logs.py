"""Check the log N log S slope of a population."""
from frbpoppy.population import unpickle

MAKE = True

if MAKE:
    from frbpoppy import CosmicPopulation, Survey, SurveyPopulation

    # Generate an FRB population
    days = 2
    population = CosmicPopulation(days*5000, days=days)

    # Setup a survey
    survey = Survey('APERTIF', gain_pattern='apertif')

    # Observe the FRB population
    surv_pop = SurveyPopulation(population, survey)
    surv_pop.save()

else:
    surv_pop = unpickle('apertif')

alpha, alpha_err, norm = surv_pop.calc_logn_logs()
print(alpha)
