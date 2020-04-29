"""Example of simulating a perfect survey."""
import matplotlib.pyplot as plt
import numpy as np
from frbpoppy import CosmicPopulation, Survey, SurveyPopulation
from convenience import hist, plot_aa_style


# Generate an FRB population
cosmic_pop = CosmicPopulation.simple(1e4, generate=True)
cosmic_pop.set_lum(model='powerlaw', low=1e40, high=1e45)

# Setup a survey
survey = Survey('perfect')
survey.snr_limit = 1e9

# Observe the FRB population
survey_pop = SurveyPopulation(cosmic_pop, survey)


plot_aa_style(cols=2)

f, axes = plt.subplots(1, 2)
s = survey_pop.frbs

snr, bins = hist(s.snr, bin_type='log')
axes[0].step(np.cumsum(snr), bins, where='mid')
axes[0].set_title('SNR')

fluence, bins = hist(s.fluence, bin_type='log')
axes[1].step(np.cumsum(fluence), bins, where='mid')
axes[1].set_title(r'Fluence')

# Set axis
axes[0].set_xscale('log')
axes[0].set_yscale('log')
axes[1].set_xscale('log')
axes[1].set_yscale('log')

plt.tight_layout()
plt.show()
