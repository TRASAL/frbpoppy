"""Plot intensity profile of theoretical beam patterns."""
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic as bstat

from frbpoppy.survey import Survey

OBSERVATORIES = [('parkes', 'HTRU'),
                 ('apertif', 'APERTIF')]

n = int(1e6)

for obs in OBSERVATORIES:

    survey = obs[1]
    pattern = obs[0]

    s = Survey(survey, gain_pattern=pattern)

    data = [s.intensity_profile(test=True) for e in range(n)]

    offset, int_pro = zip(*sorted(data))

    # Offset in degrees
    offset = [o/60. for o in offset]

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
