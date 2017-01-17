import populate
import dosurvey

pop = populate.generate(1000)

survey_pop = dosurvey.run(pop, ['PMSURV'])
