# -*- coding: future_fstrings -*-
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
