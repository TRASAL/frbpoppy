# -*- coding: future_fstrings -*-
"""Calculate rates for a local population."""
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

from frbpoppy import Survey

from quick import get_cosmic_pop, get_survey_pop

MAKE = True
OBSERVE = True
ALPHAS = np.around(np.linspace(-0.2, -2.5, 7), decimals=2)
SURVEYS = ('palfa', 'htru', 'askap-fly')
SIZE = 'small'

def simple_rates(make=MAKE, observe=OBSERVE, alphas=ALPHAS, size=SIZE,
                surveys=SURVEYS, output=False):
    # Be a bit smart about which populations need to be loaded
    load = True
    if not make and not observe:
        load = False

    # Construct populations
    pops = []
    for alpha in alphas:
        pop = get_cosmic_pop('alpha_simple',
                             size,
                             load=load,
                             overwrite=make,
                             alpha=alpha)

        pops.append(pop)

    # Survey populations
    rates = defaultdict(list)

    for s in surveys:

        survey = Survey(name=s, gain_pattern='perfect', n_sidelobes=0.5)

        for i, pop in enumerate(pops):
            surv_rates = get_survey_pop(pop, survey).rates()

            if output:
                alpha = alphas[i]
                print(f'Alpha:{alpha:.2}, Survey: {s}, Det: {surv_rates.det}')

            rate = (surv_rates.det / surv_rates.days)
            rates[s].append(rate)

    # Scale rates to HTRU
    for surv in surveys:
        if surv != 'htru':
            norm = []
            for i, r in enumerate(rates[surv]):
                norm.append(r/rates['htru'][i])
            rates[surv] = norm
    rates['htru'] = [r/r for r in rates['htru']]

    return rates


def main():

    rates = simple_rates(output=True)

    for surv in rates:
        rate = rates[surv]
        plt.plot(ALPHAS, rate, label=surv)

        plt.xlabel(r'$\alpha$')
        plt.ylabel(r'Events / htru')
        plt.yscale('log')
        plt.xlim((min(ALPHAS), max(ALPHAS)))
        plt.legend()
        plt.gca().invert_xaxis()
        plt.tight_layout()
        plt.grid()
        plt.savefig('./plots/simple_rates.pdf')


if __name__ == '__main__':
    main()
