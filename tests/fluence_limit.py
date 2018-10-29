"""Test fluence limit calculation for various surveys."""
from frbpoppy import Survey

SURVEYS = ['APERTIF',
           'ASKAP-FLY',
           'CHIME',
           'FAST',
           'HTRU',
           'LOFAR',
           'MWA',
           'UTMOST']

r = '{:10.19} {:>10}'
print(r.format('Survey', 'Jy*ms'))
n = '{:10.19} {:>10.2f}'

for survey in SURVEYS:
    surv = Survey(survey)
    fluence_limit = surv.calc_fluence_limit()
    print(n.format(survey, fluence_limit))
