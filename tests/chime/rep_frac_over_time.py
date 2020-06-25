"""Input versus output poisson rate."""
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from frbpoppy import CosmicPopulation, lognormal, Survey, SurveyPopulation
from frbpoppy import unpickle

from tests.convenience import plot_aa_style, rel_path

init_surv_time_frac = 0.05
MAKE = False

N_SRCS = 3.6e4
N_DAYS = 100
RATE = 9  # per day
# Chime started in Aug 2018. Assuming 2/day for one-offs.
# Total of 9 repeaters published on 9 Aug 2019. = ~year
N_CHIME = {'rep': 9, 'one-offs': 365*2, 'time': 365}

# Set up survey
chime = Survey('chime', n_days=N_DAYS)
chime.set_beam(model='chime')

if MAKE:

    r = CosmicPopulation(N_SRCS, n_days=N_DAYS, repeaters=True)
    r.set_dist(model='vol_co', z_max=1.0)
    r.set_dm_host(model='gauss', mean=100, std=200)
    r.set_dm_igm(model='ioka', slope=1000, std=None)
    r.set_dm(mw=True, igm=True, host=True)
    r.set_emission_range(low=100e6, high=10e9)
    r.set_lum(model='powerlaw', per_source='different', low=1e40, high=1e45,
              power=0)
    r.set_si(model='gauss', mean=-1.4, std=1)
    r.set_w(model='lognormal', per_source='different', mean=0.1, std=1)
    rate = lognormal(RATE, 1, int(N_SRCS))
    r.set_time(model='poisson', rate=rate)

    # Only generate FRBs in CHIME's survey region
    r.set_direction(model='uniform',
                    min_ra=chime.ra_min,
                    max_ra=chime.ra_max,
                    min_dec=chime.dec_min,
                    max_dec=chime.dec_max)

    r.generate()

    surv_pop = SurveyPopulation(r, chime)
    surv_pop.save()
else:
    surv_pop = unpickle('cosmic_chime')

# Set up plot style
plot_aa_style(cols=1)
f, ax1 = plt.subplots(1, 1)

# See how fraction changes over time
n_pointings = chime.n_pointings
days = np.linspace(0, N_DAYS, (N_DAYS*n_pointings)+1)
fracs = []

for day in tqdm(days):
    t = surv_pop.frbs.time.copy()
    time = np.where(t < day, t, np.nan)
    n_rep = ((~np.isnan(time)).sum(1) > 1).sum()
    n_one_offs = ((~np.isnan(time)).sum(1) == 1).sum()
    frac = n_rep / (n_rep + n_one_offs)
    fracs.append(frac)

ax1.plot(days, fracs)

# Further plot details
ax1.set_xlabel(r'Time (days)')
ax1.set_ylabel(r'$N_{\textrm{repeaters}}/N_{\textrm{detections}}$')
ax1.set_xlim(0, max(days))
ax1.set_yscale('log')

# Save figure
plt.tight_layout()
plt.savefig(rel_path(f'plots/rep_frac_over_time.pdf'))
plt.clf()
