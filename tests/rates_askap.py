"""Reproduce event rates versus alpha plot from Connor et al (2017)."""
import numpy as np
from collections import defaultdict

from frbpoppy import Survey

from quick import get_cosmic_pop, get_survey_pop
from rates_alpha import plot_event_rates

MAKE = False
OBSERVE = True
PLOT = True
SIZE = 'medium'
SURVEYS = ('htru', 'askap-fly', 'askap-incoh')
ALPHAS = np.around(np.linspace(-0.2, -2.5, 7), decimals=2)


def main():
    """Run main part of code."""
    # Be a bit smart about which populations need to be loaded
    load = True
    if not MAKE and not OBSERVE:
        load = False

    # Construct populations
    pops = []
    for alpha in ALPHAS:
        pop = get_cosmic_pop('alpha',
                             SIZE,
                             load=load,
                             overwrite=MAKE,
                             alpha=alpha)
        pops.append(pop)

    # Survey populations
    plot_data = defaultdict(list)
    plot_data['alpha'] = [a for a in ALPHAS]
    for i, pop in enumerate(pops):

        alpha = ALPHAS[i]

        for s in SURVEYS:

            pattern = 'perfect'
            n_sl = 0.5

            survey = Survey(name=s, gain_pattern=pattern, n_sidelobes=n_sl)
            surv_rates = get_survey_pop(pop, survey).rates()
            print(f'Alpha:{alpha:.2}, Survey: {s}, Det: {surv_rates.det}')
            events_per_week = (surv_rates.det / surv_rates.days) * 7
            plot_data[s].append(events_per_week)

    # Plot population event rates
    if PLOT:
        plot_event_rates(plot_data, filename='plots/rates_askap.pdf',
                         plot_exp=False)


if __name__ == '__main__':
    main()
