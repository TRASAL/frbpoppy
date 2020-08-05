"""Check the log N log F slope for future surveys."""
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, hist
from frbpoppy import unpickle, pprint

from tests.convenience import plot_aa_style, rel_path

MAKE = True
SURVEYS = ('htru', 'fast', 'puma-full', 'chord', 'ska1-low', 'ska1-mid')


if MAKE:
    surv_pops = []
    pop = CosmicPopulation.complex(1e5, generate=False)
    pop.generate()

    for name in SURVEYS:
        survey = Survey(name)
        surv_pop = SurveyPopulation(pop, survey)
        surv_pop.save()
        surv_pops.append(surv_pop)
else:
    surv_pops = []
    for name in SURVEYS:
        surv_pops.append(unpickle(f'complex_{name}'))

# Start plot
plot_aa_style()
fig, ax1 = plt.subplots(1, 1)

# Fluence plot
ax1.set_xlabel('S/N')
ax1.set_xscale('log')
ax1.set_ylabel(r'\#(${>}\text{S/N}$)')
ax1.set_yscale('log')

# Update fluence plot
for i, surv_pop in enumerate(surv_pops):
    name = surv_pop.name.split('_')[-1]
    snr = surv_pop.frbs.snr

    if snr.size == 0:
        pprint(f'No FRBs in {name} population')
        continue

    bins, values = hist(snr, bin_type='log', norm=None)

    # Cumulative sum
    values = np.cumsum(values[::-1])[::-1]

    # Normalise to area on sky
    if not np.isnan(values.all()):
        values = values * surv_pop.source_rate.f_area

    plt.step(bins, values, where='mid', label=name)

plt.legend()
plt.tight_layout()
plt.savefig(rel_path('./plots/logn_logs_future_surveys.pdf'))
