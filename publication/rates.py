"""Plot rates per survey."""
from collections import namedtuple
import matplotlib.pyplot as plt
import numpy as np

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, unpickle

MAKE = False
SurveyProps = namedtuple('Survey', ['name', 'prior', 'pattern'])
surveys = []
surveys.append(SurveyProps('HTRU', 14, 'airy'))
surveys.append(SurveyProps('APERTIF', 7, 'airy'))
surveys.append(SurveyProps('UTMOST-1D', 30, 'airy'))
surveys.append(SurveyProps('ASKAP-FLY', 14, 'airy'))

if MAKE:
    n_per_day = 5000
    days = 28
    pop_std = CosmicPopulation(n_per_day*days, days=days, name='standard')
else:
    pop_std = unpickle('standard')

# Set up data for plotting
y_pos = np.arange(len(surveys))
prior = [s.prior for s in surveys]
prediction = []

for survey in surveys:
    surv = Survey(survey.name, gain_pattern=survey.pattern)
    surv_pop = SurveyPopulation(pop_std, surv)
    prediction.append(surv_pop.rates().exp)

# Create plot
fig, ax = plt.subplots()
plt.hlines(y_pos, prior, prediction, alpha=0.3, colors='k')
plt.plot(prior, y_pos, 'o', label='Prior')
plt.plot(prediction, y_pos, 'o', label='Frbpoppy')

# Add plot details
ax.set_yticks(y_pos)
ax.set_yticklabels([s.name for s in surveys])
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('Days per FRB')
plt.legend()
plt.tight_layout()
plt.xscale('log')

plt.savefig('plots/rates.pdf')
