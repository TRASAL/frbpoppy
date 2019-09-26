"""Test perfect survey."""
import matplotlib.pyplot as plt
import os

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation
from convenience import hist

# Change working directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))

cosmic_pop = CosmicPopulation.simple(1e4, generate=False)
cosmic_pop.z_max = 0.01
cosmic_pop.lum_min = 1e38
cosmic_pop.generate()

survey = Survey('perfect')

survey_pop = SurveyPopulation(cosmic_pop, survey)

# Use A&A styling for plots
plt.style.use('./aa.mplstyle')
plt.rcParams["figure.figsize"] = (5.75373, 5.75373)

f, axes = plt.subplots(2, 2)
s = survey_pop.frbs
axes[0, 0].step(*hist(s.snr, bin_type='log'), where='mid')
axes[0, 0].set_title('SNR')

axes[1, 0].step(*hist(s.T_sys, bin_type='log'), where='mid')
axes[1, 0].set_title(r'T$_{\text{sys}}$')
axes[1, 0].set_xlim(1, 1e2)
axes[1, 0].set_ylim(1e-3, 1)

axes[0, 1].step(*hist(s.s_peak, bin_type='log'), where='mid')
axes[0, 1].set_title(r'S$_{\text{peak}}$')
axes[0, 1].set_xlim(1e-3, 1)
axes[0, 1].set_ylim(1e-3, 1)

axes[1, 1].step(*hist(s.w_arr, bin_type='log'), where='mid')
axes[1, 1].set_title(r'w$_{\text{arr}}$')
axes[1, 1].set_ylim(1e-3, 1)

for x in [0, 1]:
    for y in [0, 1]:
        axes[x, y].set_xscale('log')
        axes[x, y].set_yscale('log')

plt.tight_layout()
plt.savefig('plots/perfect.pdf')
