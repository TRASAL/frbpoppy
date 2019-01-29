"""Plot intensity profile of theoretical beam patterns."""
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binned_statistic as bstat

from frbpoppy.survey import Survey

OBSERVATORIES = [('parkes', 'htru'),
                 ('apertif', 'apertif')]

n = int(1e6)

for obs in OBSERVATORIES:

    survey = obs[1]
    pattern = obs[0]

    s = Survey(survey, gain_pattern=pattern)
    int_pro, offset = s.intensity_profile(n_gen=n)

    # Sort the values
    sorted_int = np.argsort(offset)
    int_pro = int_pro[sorted_int]
    offset = offset[sorted_int]

    # Offset in degrees
    offset = offset/60.

    bins = 1e2

    bin_means, bin_edges, bin_numbers = bstat(offset,
                                              int_pro,
                                              statistic='mean',
                                              bins=bins)

    bin_mins, _, _ = bstat(offset, int_pro, statistic='min', bins=bins)
    bin_maxs, _, _ = bstat(offset, int_pro, statistic='max', bins=bins)

    center = (bin_edges[:-1] + bin_edges[1:]) / 2

    plt.plot(center, bin_means, label=pattern)
    plt.fill_between(center, bin_mins, bin_maxs, alpha=0.2)


plt.xlabel(f'Offset ($\degree$)')
plt.ylabel('Intensity Profile')
plt.yscale('log')
plt.legend()
plt.tight_layout()
plt.savefig('plots/int_pro_surveys.pdf')
