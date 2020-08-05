"""Calculate rates for a local population."""
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

from frbpoppy import (CosmicPopulation, Survey, LargePopulation, pprint,
                      unpickle)

from tests.convenience import plot_aa_style, rel_path

REMAKE = True
ALPHAS = np.around(np.linspace(-0.2, -2.5, 7), decimals=2)
SURVEYS = ('palfa', 'htru', 'askap-fly')
SIZE = 1e7


def simple_rates(remake=REMAKE, alphas=ALPHAS, size=SIZE, surveys=SURVEYS):
    """Calculate expected rates for a simple populations."""
    rates = defaultdict(list)

    # Don't always regenerate a population
    if remake is False:
        for alpha in alphas:
            for s in surveys:
                surv_rates = unpickle(f'simple_alpha_{alpha}_{s}').source_rate
                pprint(f'Alpha:{alpha:.2}, Survey: {s}, Det: {surv_rates.det}')
                rate = (surv_rates.det / surv_rates.days)
                rates[s].append(rate)
    else:
        pops = []
        for alpha in alphas:
            pop = CosmicPopulation.simple(size)
            pop.set_dist(alpha=alpha)
            pop.name = f'simple_alpha_{alpha}'
            pops.append(pop)

            # Set up surveys
            ss = []
            for s in surveys:
                survey = Survey(name=s)
                survey.set_beam(model='perfect',
                                n_sidelobes=0.5)
                ss.append(survey)

            surv_pops = LargePopulation(pop, *ss).pops

            for i, s in enumerate(surveys):
                surv_rates = surv_pops[i].source_rate
                pprint(f'Alpha:{alpha:.2}, Survey: {s}, Det: {surv_rates.det}')
                rate = (surv_rates.det / surv_rates.days)
                rates[s].append(rate)

    # Scale rates to HTRU
    for s in surveys:
        if s != 'htru':
            norm = []
            for i, r in enumerate(rates[s]):
                norm.append(r/rates['htru'][i])
            rates[s] = norm
    rates['htru'] = [r/r for r in rates['htru']]

    return rates


def main():
    """Plot expected simple rates."""
    plot_aa_style()

    rates = simple_rates()
    for surv in rates:
        rate = rates[surv]
        plt.plot(ALPHAS, rate, label=surv)

    plt.xlabel(r'$\alpha_{\text{in}}$')
    plt.ylabel(r'Events / htru')
    plt.yscale('log')
    plt.xlim((min(ALPHAS), max(ALPHAS)))
    plt.legend()
    plt.gca().invert_xaxis()
    plt.tight_layout()
    plt.grid()
    plt.savefig(rel_path('./plots/rates_simple.pdf'))


if __name__ == '__main__':
    main()
