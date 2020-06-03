"""Check the log N log F slope of a population."""
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, hist
from frbpoppy import unpickle

from tests.convenience import plot_aa_style, rel_path

MAKE = True
SURVEYS = ('htru', 'fast', 'puma', 'chord')


if MAKE:
    surv_pops = []
    pop = CosmicPopulation.complex(1e8, generate=True)
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
# ax1.set_xlim(9, 1e5)
# ax1.set_ylim(1e-3, 1e3)

# Update fluence plot
for i, surv_pop in enumerate(surv_pops):
    name = surv_pop.name.split('_')[-1]
    snr = surv_pop.frbs.snr
    try:
        bins, values = hist(snr, bin_type='log', norm=None)
    except ValueError:
        bins, values = np.array([np.nan]), np.array([np.nan])

    # Cumulative sum
    values = np.cumsum(values[::-1])[::-1]

    # Normalise to area on sky
    if not np.isnan(values.all()):
        values = values * surv_pop.source_rate.f_area

    plt.step(bins, values, where='mid', label=name)

plt.legend()
plt.tight_layout()
plt.savefig(rel_path('./plots/logn_logs_future_surveys.pdf'))
