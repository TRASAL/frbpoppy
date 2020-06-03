"""Show how logN-LogS differs per survey."""
import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from frbpoppy import paths, unpickle, hist
from tests.convenience import plot_aa_style, rel_path

SURVEYS = ('askap-fly', 'fast', 'htru', 'apertif', 'puma', 'chord')
vars = {'alpha': -1.5, 'li': 0.0, 'si': 0.0}


def get_pops(alpha='*', li='*', si='*', survey='*'):
    filename = f'complex_alpha_{alpha}_lum_{li}_si_{si}_{survey}.p'
    filter = os.path.join(paths.populations(), filename)
    pop_paths = glob(filter)

    pops = []
    for path in pop_paths:
        if '_for_plotting' not in path:
            pops.append(unpickle(path))
    return pops


# Start plot
plot_aa_style()
fig, ax1 = plt.subplots(1, 1)

# Fluence plot
ax1.set_xlabel('S/N')
ax1.set_xscale('log')
ax1.set_ylabel(r'N $>$ S/N')
ax1.set_yscale('log')
ax1.set_xlim(9, 1e5)
ax1.set_ylim(1e-3, 1e3)

# Update fluence plot
for i, survey in enumerate(SURVEYS):
    pop = get_pops(**vars, survey=survey)[0]

    snr = pop.frbs.snr
    try:
        bins, values = hist(snr, bin_type='log', norm=None)
    except ValueError:
        bins, values = np.array([np.nan]), np.array([np.nan])

    # Cumulative sum
    values = np.cumsum(values[::-1])[::-1]
    # Normalise to area on sky
    values = values * pop.source_rate.f_area

    plt.step(bins, values, where='mid', label=survey)

plt.legend()
plt.tight_layout()
plt.savefig(rel_path('./plots/logn_logs_surveys.pdf'))
