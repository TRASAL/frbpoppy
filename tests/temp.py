"""Survey brightness distribution."""
from frbpoppy import Survey, SurveyPopulation, plot, unpickle

# pop = unpickle('alpha_large_-2.5')
# surv = Survey('perfect', gain_pattern='perfect')
# surv_pop = SurveyPopulation(pop, surv)
# surv_pop.name = 'alpha_large_-2.5_perfect'
# surv_pop.save()

# perfect = unpickle('alpha_large_-2.5_perfect')
# palfa = unpickle('alpha_large_-2.5_palfa')
# plot(perfect, palfa, frbcat=False)

large = unpickle('alpha_large_-2.5_palfa')

plot(large, frbcat=False)
