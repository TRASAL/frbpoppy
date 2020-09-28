"""Check the log N log F slope of various populations."""
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, hist
from frbpoppy import unpickle

from tests.convenience import plot_aa_style, rel_path

MAKE = True
SURVEYS = ('askap-fly', 'fast-crafts', 'parkes-htru', 'wrst-apertif')


if MAKE:
    surv_pops = []
    pop = CosmicPopulation.simple(1e7, generate=True)
    pop.set_dist(model='vol_co', z_max=0.01)
    pop.set_lum(model='constant', value=1e36)
    pop.set_dm(mw=False, igm=True, host=True)
    pop.set_w(model='lognormal', mean=0.1, std=0.7)
    pop.name = 'local'
    pop.generate()

    for name in SURVEYS:
        survey = Survey(name)
        surv_pop = SurveyPopulation(pop, survey)
        surv_pop.save()
        surv_pops.append(surv_pop)
else:
    surv_pops = []
    for name in SURVEYS:
        surv_pops.append(unpickle(f'local_{name}'))

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
plt.savefig(rel_path('./plots/logn_logs_local_lum_1e36_w_int.pdf'))
