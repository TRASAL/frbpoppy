"""Plot rates per survey."""
from collections import namedtuple
import matplotlib.pyplot as plt
import numpy as np

from frbpoppy.do_populate import generate
from frbpoppy.do_survey import observe
from frbpoppy.population import unpickle

MAKE = False
Survey = namedtuple('Survey', ['name', 'prior', 'pattern'])
surveys = []
surveys.append(Survey('HTRU', 14, 'airy'))
surveys.append(Survey('APERTIF', 7, 'airy'))
surveys.append(Survey('UTMOST-1D', 30, 'airy'))
surveys.append(Survey('ASKAP-FLY', 14, 'airy'))

if MAKE:
    n_per_day = 5000
    days = 28
    pop_std = generate(n_per_day*days, days=days, name='standard')
else:
    pop_std = unpickle('standard')

# Set up data for plotting
y_pos = np.arange(len(surveys))
prior = [s.prior for s in surveys]
prediction = []

for survey in surveys:
    survey = observe(pop_std,
                     survey.name,
                     gain_pattern=survey.pattern,
                     return_survey=True,
                     return_pop=False,
                     output=True)

    prediction.append(survey.days_per_frb)

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
# plt.xscale('log')

plt.savefig('plots/rates.pdf')
