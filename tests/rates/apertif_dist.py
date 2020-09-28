"""Calculate the expected detection rates for apertif."""
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, hist

from tests.convenience import plot_aa_style, rel_path
from alpha_real import EXPECTED, poisson_interval

N_DAYS = 1  # Not used in eventual result
SCALE_TO = 'parkes-htru'

pop = CosmicPopulation.complex(n_srcs=1e5, n_days=N_DAYS)
pop.generate()

apertif = Survey('wrst-apertif', n_days=N_DAYS)
apertif.set_beam(model='apertif_real')

if SCALE_TO == 'parkes-htru':
    htru = Survey('parkes-htru', n_days=N_DAYS)
    htru.set_beam(model='parkes')
if SCALE_TO == 'askap':
    askap = Survey('askap-fly', n_days=N_DAYS)
    askap.set_beam(model='gaussian', n_sidelobes=0.5)

days_per_frbs = []
for i in tqdm(range(2000), desc='Survey Run'):

    apertif_pop = SurveyPopulation(pop, apertif, mute=True)

    if SCALE_TO == 'parkes-htru':
        htru_pop = SurveyPopulation(pop, htru, mute=True)
        n_frbs_htru = EXPECTED['parkes-htru'][0]
        n_days_htru = EXPECTED['parkes-htru'][1]
        scaled_n_days = n_days_htru*(htru_pop.source_rate.det / n_frbs_htru)

    if SCALE_TO == 'askap':
        askap_pop = SurveyPopulation(pop, askap, mute=True)
        n_frbs_askap = EXPECTED['askap-fly'][0]
        n_days_askap = EXPECTED['askap-fly'][1]
        scaled_n_days = n_days_askap*(askap_pop.source_rate.det / n_frbs_askap)

    days_per_frb = scaled_n_days / apertif_pop.source_rate.det
    # print(f'{days_per_frb} days per frb')
    days_per_frbs.append(days_per_frb)


days_per_frbs = np.array(days_per_frbs)
mean = np.mean(days_per_frbs)
std = np.std(days_per_frbs)
poisson_std = poisson_interval(mean)
print(f'Mean rate for apertif is {mean}')
print(f'Standard deviation of {std}')

# Plot
rates, values = hist(days_per_frbs, bin_type='lin')
plot_aa_style()
plt.step(rates, values, where='mid')
plt.xlabel(f'Apertif days per burst scaled to {SCALE_TO}')
plt.savefig(rel_path('./plots/rate_apertif_dist.pdf'))
