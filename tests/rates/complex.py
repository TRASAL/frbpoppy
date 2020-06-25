"""Try creating the most realistic event rates."""
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

from frbpoppy import (CosmicPopulation, Survey, LargePopulation, pprint,
                      unpickle)

from tests.convenience import plot_aa_style, rel_path

REMAKE = True
ALPHAS = np.around(np.linspace(-0.5, -2.0, 7), decimals=2)
SURVEYS = ('askap-fly', 'fast', 'htru', 'apertif', 'palfa')
SIZE = 1e8


def complex_rates(remake=REMAKE, alphas=ALPHAS, size=SIZE, surveys=SURVEYS):
    """Calculate expected rates for a complex populations."""
    rates = defaultdict(list)

    # Don't always regenerate a population
    if remake is False:
        for alpha in alphas:
            for s in surveys:
                surv_rates = unpickle(f'complex_alpha_{alpha}_{s}').source_rate
                pprint(f'Alpha:{alpha:.2}, Survey: {s}, Det: {surv_rates.det}')
                rate = (surv_rates.det / surv_rates.days)
                rates[s].append(rate)
    else:
        pops = []
        for alpha in alphas:
            if alpha <= -1.0:
                size = 1e7
            if alpha <= -1.5:
                size = 1e8
            pop = CosmicPopulation.complex(size)
            pop.set_dist(model='vol_co', z_max=2.5, alpha=alpha,
                         H_0=67.74, W_m=0.3089, W_v=0.6911)
            pop.set_lum(model='powerlaw', low=1e40, high=1e45, power=-1)
            pop.name = f'complex_alpha_{alpha}'
            pops.append(pop)

            # Set up surveys
            ss = []
            for s in surveys:
                survey = Survey(name=s)
                survey.set_beam(model='airy', n_sidelobes=1)
                ss.append(survey)

            surv_pops = LargePopulation(pop, *ss).pops

            for i, s in enumerate(surveys):
                surv_rates = surv_pops[i].source_rate
                pprint(f'Alpha:{alpha:.2}, Survey: {s}, Det: {surv_rates.det}')
                rate = (surv_rates.det / surv_rates.days)
                rates[s].append(rate)

    # Scale rates to first survey in list
    for s in surveys:
        if s != surveys[0]:
            norm = []
            for i, r in enumerate(rates[s]):
                norm.append(r/rates[surveys[0]][i])
            rates[s] = norm
    rates[surveys[0]] = [r/r for r in rates[surveys[0]]]

    return rates


def main():
    """Plot expected complex rates."""
    plot_aa_style()

    rates = complex_rates()
    for surv in rates:
        rate = rates[surv]
        plt.plot(ALPHAS, rate, label=surv)

    plt.xlabel(r'$\alpha_{\text{in}}$')
    plt.ylabel(rf'Events / {SURVEYS[0]}')
    plt.yscale('log')
    plt.xlim((min(ALPHAS), max(ALPHAS)))
    plt.legend()
    plt.gca().invert_xaxis()
    plt.tight_layout()
    plt.grid()
    plt.savefig(rel_path('./plots/complex_rates.pdf'))


if __name__ == '__main__':
    main()
