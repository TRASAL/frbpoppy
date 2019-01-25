"""Reproduce event rates versus alpha plot from Connor et al (2017)."""
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pickle

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, unpickle

SURVEYS = ('htru', 'askap-fly', 'askap-incoh')
ALPHAS = np.around(np.linspace(-0.2, -2.5, 7), decimals=2)
MAKE = False
OBSERVE = False

plot_data = defaultdict(list)
pops = []


def save_obj(obj, name):
    with open('obj/' + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


if MAKE:
    for alpha in ALPHAS:
        n_per_day = 5000
        days = 200

        # Local, standard candle population with constant pulse widths
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
                               si_mu=0.,
                               si_sigma=0.,
                               z_max=0.1)

        pop.save()
        pops.append(pop)
        plot_data['alpha'].append(alpha)

if OBSERVE:

    for i, alpha in enumerate(ALPHAS):

        if not MAKE:
            pop = unpickle(f'alpha-{alpha}')
        else:
            pop = pops[i]

        for s in SURVEYS:

            pattern = 'perfect'
            survey = Survey(name=s, gain_pattern=pattern, n_sidelobes=0.5)
            surv_rates = SurveyPopulation(pop, survey).rates()
            print(f'Alpha:{pop.alpha:.2}, Survey: {s}, Det: {surv_rates.det}')
            events_per_week = (surv_rates.det / surv_rates.days) * 7
            events_per_week *= 1e3  # Just upping the values
            plot_data[s].append(events_per_week)

    save_obj(plot_data, 'rates')
else:
    plot_data = load_obj('rates')

plot_data['alpha'] = [a for a in ALPHAS]

for survey in plot_data:
    if survey == 'alpha':
        continue
    ps = np.array(plot_data[survey])
    ph = np.array(plot_data['htru'])
    pa = np.array(plot_data['alpha'])

    if survey != 'HTRU':
        ps[ph == 0.] = 0.

    # Normalise wrt to Parkes at all alphas
    norm = np.array([ps[i]/ph[i] for i in range(len(ps))])
    norm[np.isnan(norm)] = 0.
    cleaned_norm = norm[norm > 0.]
    alpha = pa[norm > 0.]

    plt.plot(alpha, cleaned_norm, marker='o', label=survey.title())

plt.xlabel(r'$\alpha$')
plt.ylabel('Events / Htru')
plt.yscale('log')
plt.gca().invert_xaxis()
plt.legend()#loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.savefig('plots/rates_askap.pdf')
