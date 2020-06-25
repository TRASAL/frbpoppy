"""Check the log N log F slope of a population."""
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import Frbcat, hist

from tests.convenience import plot_aa_style, rel_path

SURVEYS = ('askap-fly', 'crafts', 'htru', 'apertif', 'askap-incoh')

frbcat = Frbcat().df
snrs = {s: frbcat[frbcat.survey == s].snr.values for s in SURVEYS}
snrs['apertif'] = [15.4, 12.9, 13.2, 60, 13, 38.1, 18.1, 27, 16.7, 13, 18.8, 27.7, 18.9, 17.84, 10.2, 14.84, 10.25]
snrs['crafts'] = [19]
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
