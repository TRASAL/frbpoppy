"""Check the current state of the log N log F."""
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import Frbcat, hist

from tests.convenience import plot_aa_style, rel_path

SURVEYS = ('askap-fly', 'fast-crafts', 'parkes-htru', 'wrst-apertif',
           'askap-incoh')

frbcat = Frbcat().df
snrs = {s: frbcat[frbcat.survey == s].snr.values for s in SURVEYS}
snrs['wrst-apertif'] = [15.4, 12.9, 13.2, 60, 13, 38.1, 18.1, 27, 16.7, 13, 18.8,
                   27.7, 18.9, 17.84, 10.2, 14.84, 10.25]
snrs['fast-crafts'] = [19]
# Start plot
plot_aa_style()
fig, ax1 = plt.subplots(1, 1)

# Fluence plot
ax1.set_xlabel('S/N')
ax1.set_xscale('log')
ax1.set_ylabel(r'\#(${>}\text{S/N}$)')
ax1.set_yscale('log')

# Update fluence plot
for survey in SURVEYS:
    name = survey
    snr = snrs[survey]

    try:
        bins, values = hist(snr, bin_type='log', norm=None)

        # Cumulative sum
        values = np.cumsum(values[::-1])[::-1]

        plt.step(bins, values, where='mid', label=name)
    except ValueError:
        continue

plt.legend()
plt.tight_layout()
plt.savefig(rel_path('./plots/logn_logs_current.pdf'))
