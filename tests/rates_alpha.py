"""Reproduce event rates versus alpha plot from Connor et al (2017)."""
from collections import defaultdict
from scipy.stats import chi2, norm
from scipy.integrate import quad
import matplotlib.pyplot as plt
import numpy as np

from frbpoppy import Survey

from quick import get_cosmic_pop, get_survey_pop

MAKE = True
OBSERVE = True
PLOT = True
SIZE = 'medium'
SURVEYS = ('palfa', 'htru', 'askap-fly')
ALPHAS = np.around(np.linspace(-0.2, -2.5, 7), decimals=2)

expected = {'htru': [9, 24 * 0.551 / 1549],  # N_frbs, scaling to get frbs/day
            'apertif': [1, 1 / 7],
            'askap-fly': [20, 24 / 32840 * 8],
            'utmost': [0, 0],
            'chime': [0, 0],
            'palfa': [1, 1 / 24.1],
            'guppi': [0.4, 1 / 81]  # 0.4 is my own assumption
            }


def poisson_interval(k, sigma=1):
    """
    Use chi-squared info to get the poisson interval.

    Give a number of observed events, which range of observed events would have
    been just as likely given a particular interval?

    Based off https://stackoverflow.com/questions/14813530/
    poisson-confidence-interval-with-numpy
    """
    gauss = norm(0, 1).pdf
    a = 1 - quad(gauss, -sigma, sigma, limit=1000)[0]
    low, high = (chi2.ppf(a/2, 2*k) / 2, chi2.ppf(1-a/2, 2*k + 2) / 2)
    if k == 0:
        low = 0.0

    return low, high


def plot_event_rates(data, filename=None, plot_exp=True):
    """Plot event rates.

    Args:
        data (dict): Dictionary with per survey the event rates

    """
    # Use a nice font for axes
    plt.rc('text', usetex=True)

    f, ax = plt.subplots(1, 1)

    a = 0.25
    for survey in data:
        if survey == 'alpha':
            continue

        color = next(ax._get_lines.prop_cycler)['color']

        ps = np.array(data[survey])
        ph = np.array(data['htru'])
        pa = np.array(data['alpha'])

        if survey != 'htru':
            ps[ph == 0.] = 0.

        # Normalise wrt to Parkes at all alphas
        norm = np.array([ps[i]/ph[i] for i in range(len(ps))])
        norm[np.isnan(norm)] = 0.
        cleaned_norm = norm[norm > 0.]
        alpha = pa[norm > 0.]

        ax.plot(alpha, cleaned_norm, marker='o', label=survey, color=color)

        if plot_exp:

            # Plot expected rate
            exp_n = expected[survey][0]
            exp_scaling = expected[survey][1]

            norm = 1 / (expected['htru'][0] * expected['htru'][1])

            exp_min, exp_max = poisson_interval(exp_n, sigma=2)

            exp = exp_n * exp_scaling * norm
            exp_min *= exp_scaling * norm
            exp_max *= exp_scaling * norm

            ax.plot(pa, [exp]*len(pa), ls='--', color=color, alpha=a+0.2)
            ax.fill_between(pa, [exp_min]*len(pa), [exp_max]*len(pa),
                            alpha=a, color=color)
            a += 0.1

    plt.xlabel(r'$\alpha$')
    plt.ylabel('Events / Htru')
    plt.yscale('log', nonposy='mask')
    plt.gca().invert_xaxis()
    plt.legend()
    plt.tight_layout()

    if filename is None:
        filename = 'plots/event_rates_boxcar.pdf'

    plt.savefig(filename)


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
        plot_event_rates(plot_data)


if __name__ == '__main__':
    main()
