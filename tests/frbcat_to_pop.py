"""Convert a Pandas DataFrame with Frbcat to a Population class."""
from frbpoppy.frbcat import Frbcat
from frbpoppy import Survey, SurveyPopulation

frbcat = Frbcat().to_pop()
survey = Survey('APERTIF')
surv_pop = SurveyPopulation(frbcat, survey)

names = surv_pop.get('name')
fluences = surv_pop.get('fluence')

for f in fluences:
    print(f)
