"""Compare the resulting logNlogS and logNLog(S/N)."""
import matplotlib.pyplot as plt
import numpy as np
from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, hist
from tests.convenience import plot_aa_style


# Generate an FRB population
cosmic_pop = CosmicPopulation.simple(1e4)
cosmic_pop.set_lum(model='powerlaw', low=1e40, high=1e45)
cosmic_pop.generate()

# Setup a survey
survey = Survey('perfect')

# Observe the FRB population
survey_pop = SurveyPopulation(cosmic_pop, survey)


plot_aa_style(cols=2)

f, axes = plt.subplots(1, 2)
s = survey_pop.frbs

bins, values = hist(s.snr, bin_type='log')
# Cumulative sum
values = np.cumsum(values[::-1])[::-1]
axes[0].step(bins, values, where='mid')
axes[0].set_title('SNR')

bins, values = hist(s.fluence, bin_type='log')
# Cumulative sum
values = np.cumsum(values[::-1])[::-1]
axes[1].step(bins, values, where='mid')
axes[1].set_title(r'Fluence')

# Set axis
axes[0].set_xscale('log')
axes[0].set_yscale('log')
axes[1].set_xscale('log')
axes[1].set_yscale('log')

plt.tight_layout()
plt.show()
