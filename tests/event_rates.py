"""Reproduce detection rates as shown in Fig. 1 from Connor (2017)."""
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, unpickle

SURVEYS = ('HTRU', 'APERTIF', 'UTMOST', 'ASKAP-FLY')
ALPHAS = (1.5001, -1.0, -1.5)
SPECTRAL_INDEX = 0.
MAKE = True
OBSERVE = True

k = 0

if MAKE:

    for alpha in ALPHAS:

        # For testing
        if k >= 1:
            continue
        k += 1

        n_per_day = 5000
        days = 3

        pop = CosmicPopulation(n_per_day*days,
                               days=days,
                               name=f'alpha-{alpha}',
                               H_0=69.6,
                               W_m=0.286,
                               W_v=0.714,
                               dm_host_model='normal',
                               dm_host_mu=100.,
                               dm_host_sigma=0.,
                               dm_igm_index=1200,
                               dm_igm_sigma=0.,
                               dm_mw_model='ne2001',
                               emission_range=[10e6, 10e9],
                               lum_range=[1e40, 1e40],
                               lum_index=0,
                               n_model='vol_co',
                               alpha=alpha,
                               pulse_model='uniform',
                               pulse_range=[1., 1.],
                               pulse_mu=1.,
                               pulse_sigma=0.,
                               repeat=0.,
                               si_mu=SPECTRAL_INDEX,
                               si_sigma=0.,
                               z_max=0.1)
        pop.save()
else:
    pop = unpickle('alpha-1.0')

if OBSERVE:

    survey = Survey('PERFECT', gain_pattern='perfect', sidelobes=0)
    surv_pop = SurveyPopulation(pop, survey)
    surv_pop.name = 'testing-alpha'
    surv_pop.save()
else:
    surv_pop = unpickle('testing-alpha')

fluences = np.array(surv_pop.get('fluence'))

# N(S)
number, bins = np.histogram(np.log10(fluences), bins=500)
# N(>S) from N(S)
n_gt_s = np.cumsum(number[::-1])[::-1]
# logS
x = bins[:-1]
# log(N(>S))
y = np.log10(n_gt_s)
# Calculate derivative
der = np.diff(y) / np.diff(x)
bin_centres = (x[:-1] + x[1:]) / 2

plt.plot(bin_centres, der)
plt.xlabel('log S')
plt.ylabel(r'$\alpha$')
plt.ylim(np.mean(der)-2, np.mean(der)+2)

plt.tight_layout()
plt.savefig('plots/event_rates.pdf')

plt.clf()


plt.plot(x, y)
plt.xlabel('log S')
plt.ylabel(r'log N(>S)')
plt.tight_layout()
plt.savefig('plots/logNlogS.pdf')
