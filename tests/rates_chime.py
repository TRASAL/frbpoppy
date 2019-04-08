"""Event rates for CHIME versus spectral index over alpha."""
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

from frbpoppy import Survey

from quick import get_cosmic_pop, get_survey_pop

MAKE = False
OBSERVE = True
PLOT = True
SIZE = 'small'
SURVEY = 'chime'
ALPHAS = np.around(np.linspace(-0.2, -2.5, 7), decimals=2)
GAMMAS = [-2., 0, 2]


def plot_rates(rates, filename='plots/rates_chime.pdf'):
    """Plot the event rates.

    Args:
        rates (dict): Dictionary of rates, with every spectral index key
            mapping to a list of event rates per alpha.
        filename (str): Where to save the file

    """
    for gamma in rates:

        rate = np.array(rates[gamma])
        xs = ALPHAS[~np.isinf(np.log10(rate))]
        rate = rate[~np.isinf(np.log10(rate))]

        plt.plot(xs, rate, marker='o', label=fr'$\gamma =$ {int(gamma)}')

    plt.xlabel(r'$\alpha$')
    plt.ylabel('Event rate')
    plt.yscale('log')
    plt.gca().invert_xaxis()
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename)


def main():
    """Run main part of code."""
    # Be a bit smart about which populations need to be loaded
    load = True
    if not MAKE and not OBSERVE:
        load = False

    # Construct populations
    pops = {}
    for gamma in GAMMAS:
        pops[gamma] = []
        for alpha in ALPHAS:
            pop = get_cosmic_pop('alpha_gamma',
                                 SIZE,
                                 load=load,
                                 overwrite=MAKE,
                                 alpha=alpha,
                                 gamma=gamma)
            pops[gamma].append(pop)

    # Survey populations
    rates = defaultdict(list)
    for gamma in GAMMAS:
        print(gamma)
        for alpha in ALPHAS:

            pattern = 'perfect'
            n_sl = 0.5
            survey = Survey(name=SURVEY, gain_pattern=pattern,
                            n_sidelobes=n_sl)
            surv_rates = get_survey_pop(pop, survey).rates()

            m = f'Alpha:{alpha:.2}, Gamma:{gamma}, Det: {surv_rates.det}'
            print(m)

            events = (surv_rates.det / surv_rates.days)
            rates[gamma].append(events)

    # Plot population event rates
    if PLOT:
        plot_rates(rates, filename='plots/rates_chime.pdf')


if __name__ == '__main__':
    main()


#
# SURVEY = 'chime'
# GAMMAS = [-2., 0, 2]
# ALPHAS = np.around(np.linspace(-0.2, -2.0, 6), decimals=2)
# MAKE = False
# OBSERVE = False
#
# rates = defaultdict(list)
# pops = []
# gammas = []
# alphas = []
#
# def save_obj(obj, name):
#     with open('obj/' + name + '.pkl', 'wb') as f:
#         pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
#
#
# def load_obj(name):
#     with open('obj/' + name + '.pkl', 'rb') as f:
#         return pickle.load(f)
#
#
# for gamma in GAMMAS:
#     for alpha in ALPHAS:
#         if MAKE:
#             n_per_day = 5000
#             days = 50
#
#             # Local, standard candle population with constant pulse widths
#             pop = CosmicPopulation(n_per_day*days,
#                                    days=days,
#                                    name=f'alpha-{alpha}-gamma-{gamma}',
#                                    H_0=67.74,
#                                    W_m=0.3089,
#                                    W_v=0.6911,
#                                    dm_host_model='normal',
#                                    dm_host_mu=100.,
#                                    dm_host_sigma=0.,
#                                    dm_igm_index=1200,
#                                    dm_igm_sigma=0.,
#                                    dm_mw_model='ne2001',
#                                    emission_range=[10e6, 10e9],
#                                    lum_range=[1e40, 1e40],
#                                    lum_index=0,
#                                    n_model='vol_co',
#                                    alpha=alpha,
#                                    pulse_model='uniform',
#                                    pulse_range=[1., 1.],
#                                    pulse_mu=1.,
#                                    pulse_sigma=0.,
#                                    si_mu=gamma,
#                                    si_sigma=0.,
#                                    z_max=2.5)
#
#             pop.save()
#             pops.append(pop)
#         alphas.append(alpha)
#         gammas.append(gamma)
#
#
# if OBSERVE:
#
#     for i, alpha in enumerate(alphas):
#
#         gamma = gammas[i]
#
#         if not MAKE:
#             pop = unpickle(f'alpha-{alpha}-gamma-{gamma}')
#         else:
#             pop = pops[i]
#
#         pattern = 'perfect'
#         survey = Survey(name=SURVEY, gain_pattern=pattern, n_sidelobes=0.5)
#         surv_rates = SurveyPopulation(pop, survey).rates()
#         print(f'Alpha:{pop.alpha:.2}, Survey: {SURVEY}, Det: {surv_rates.det}')
#         events_per_week = (surv_rates.det / surv_rates.days) * 7
#         events_per_week *= 10000  # For plotting reasons
#         rates[gamma].append(events_per_week)
#
#     save_obj(rates, 'rates')
# else:
#     rates = load_obj('rates')
#
# for gamma in rates:
#
#     rate = np.array(rates[gamma])
#
#     xs = ALPHAS[~np.isinf(np.log10(rate))]
#     rate = rate[~np.isinf(np.log10(rate))]
#
#     plt.plot(xs, rate, marker='o', label=fr'$\gamma =$ {int(gamma)}')
#
# plt.xlabel(r'$\alpha$')
# plt.ylabel('Event rate')
# plt.yscale('log')
# plt.gca().invert_xaxis()
# plt.legend()
# plt.tight_layout()
# plt.savefig('plots/rates_chime.pdf')
