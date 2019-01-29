"""Reproduce event rates versus alpha plot from Connor et al (2017)."""
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pickle

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, unpickle

SURVEYS = ('htru', 'askap-fly', 'palfa')
ALPHAS = np.around(np.linspace(-0.2, -2.5, 7), decimals=2)
MAKE = False
OBSERVE = False

expected = {'htru': 9 * 24 * 0.551 / 1549, #  frbs/day
            'apertif': 1 / 7,
            'askap-fly': 20 * 24 / 32840,
            'utmost': 0,
            'chime': 0,
            'palfa': 1 / 24.1,
            'guppi': 0.4 / 81 #  0.4 is my own assumption
}

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
        n = 1000000
        days = 1

        # Local, standard candle population with constant pulse widths
        pop = CosmicPopulation(n,
                               days=days,
                               name=f'alpha-{alpha}',
                               H_0=69.6,
                               W_m=0.286,
                               W_v=0.714,
                               dm_host_model='normal',
                               dm_host_mu=100.,
                               dm_host_sigma=0.,
                               dm_igm_index=1000,
                               dm_igm_sigma=None,
                               dm_mw_model='ne2001',
                               emission_range=[10e6, 10e9],
                               lum_range=[1e40, 1e45],
                               lum_index=0,
                               n_model='vol_co',
                               alpha=alpha,
                               pulse_model='uniform',
                               pulse_range=[1., 1.],
                               pulse_mu=1.6,
                               pulse_sigma=1.,
                               si_mu=-1.4,
                               si_sigma=1.,
                               z_max=2.5)

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

            pattern = 'airy'
            n_sl = 1
            survey = Survey(name=s, gain_pattern=pattern, n_sidelobes=n_sl)
            surv_rates = SurveyPopulation(pop, survey).rates()
            print(f'Alpha:{pop.alpha:.2}, Survey: {s}, Det: {surv_rates.det}')
            events_per_week = (surv_rates.det / surv_rates.days) * 7
            plot_data[s].append(events_per_week)

    save_obj(plot_data, 'rates')
else:
    plot_data = load_obj('rates')

# Plot the results

plot_data['alpha'] = [a for a in ALPHAS]

f, ax = plt.subplots(1, 1)

for survey in plot_data:
    if survey == 'alpha':
        continue

    ps = np.array(plot_data[survey])
    ph = np.array(plot_data['htru'])
    pa = np.array(plot_data['alpha'])

    if survey != 'htru':
        ps[ph == 0.] = 0.

    # Normalise wrt to Parkes at all alphas
    norm = np.array([ps[i]/ph[i] for i in range(len(ps))])
    norm[np.isnan(norm)] = 0.
    cleaned_norm = norm[norm > 0.]
    alpha = pa[norm > 0.]

    ax.plot(alpha, cleaned_norm, marker='o', label=survey)

    # Plot expected rate
    exp = expected[survey] / expected['htru']

    exp_min = exp - np.sqrt(exp)
    exp_max = exp + np.sqrt(exp)


    if exp == 0:
        exp_min = 1
        exp_max = 1

    if exp_min <= 0:
        y_min = 1e-2
        exp_min = y_min

    if exp_min > exp_max:
        exp_min, exp_max = exp_max, exp_min

    ax.fill_between(pa, [exp_min]*len(pa), [exp_max]*len(pa), alpha=0.3)

plt.xlabel(r'$\alpha$')
plt.ylabel('Events / Htru')
plt.yscale('log', nonposy='mask')
plt.gca().invert_xaxis()
plt.legend()
plt.tight_layout()
plt.savefig('plots/event_rates.pdf')
