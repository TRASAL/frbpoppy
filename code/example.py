import populate
import dosurvey

pop = populate.generate(10000, electron_model='ne2001')
pop.write_ascii()
survey_pop = dosurvey.run(pop, ['PMSURV'])
