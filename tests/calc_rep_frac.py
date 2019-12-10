"""Calculate fraction of repeaters in detected frbs."""
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import RepeaterPopulation, Survey, SurveyPopulation

from convenience import plot_aa_style, rel_path

MAX_DAYS = 365

r = RepeaterPopulation.simple(1e4)
r.lum_min = 1e40
r.lum_max = 1e45
r.lum_pow = 0
r.lum_rep_model = 'independent'
r.z_max = 0.01
r.times_rep_model = 'even'
r.n_days = MAX_DAYS*2

# Set DM distributions
r.dm_host_model = 'gaussian'
r.dm_host_mu = 0
r.dm_host_sigma = 0
r.dm_igm_index = 1000
r.dm_igm_sigma = 0
r.dm_mw_model = 'zero'
r.generate()

survey = Survey('perfect-small', strategy='regular')
survey.n_days = MAX_DAYS
survey.beam_pattern = 'chime'
survey.snr_limit = 1e12

surv_pop = SurveyPopulation(r, survey)

# See how fraction changes over time
days = np.linspace(0, MAX_DAYS, MAX_DAYS+1)
fracs = []
for day in days:
    t = surv_pop.frbs.time.copy()
    time = np.where(t < day, t, np.nan)

    n_rep = ((~np.isnan(time)).sum(1) > 1).sum()
    n_one_offs = ((~np.isnan(time)).sum(1) == 1).sum()
    frac = n_rep / (n_rep + n_one_offs)
    fracs.append(frac)

plot_aa_style(cols=2)
plt.plot(days, fracs)
plt.xlabel(r'Time (days)')
plt.ylabel(r'$N_{\textrm{repeaters}}/N_{\textrm{detections}}$')
plt.ylim(0, 1)
plt.tight_layout()
plt.savefig(rel_path(f'plots/rep_frac.pdf'))
plt.clf()
