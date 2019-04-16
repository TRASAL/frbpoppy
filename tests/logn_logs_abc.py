"""Plot a log N / log S graph for three different populations."""
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import Survey, SurveyPopulation

from quick import get_cosmic_pop

MAKE = False
SIZE = 'small'
GAMMAS = [-1.4, 1]


def get_local(make=MAKE, size=SIZE):

    # Construct populations
    pop = get_cosmic_pop('alpha_simple',
                         size,
                         load=True,
                         overwrite=make,
                         alpha=-1.5)

    # Survey populations
    survey = Survey(name='perfect', gain_pattern='perfect', n_sidelobes=0.5)
    surv_pop = SurveyPopulation(pop, survey)

    return surv_pop.frbs.fluence


def get_further(gamma, make=MAKE, size=SIZE):
    """Construct populations going further out."""
    # Construct populations
    pop = get_cosmic_pop('gamma',
                         size,
                         load=True,
                         overwrite=make,
                         gamma=gamma)

    if gamma == 1:
        pop.frbs.lum_bol = np.ones_like(pop.frbs.lum_bol)*10**43

    # Survey populations
    survey = Survey(name='perfect', gain_pattern='perfect', n_sidelobes=0.5)
    surv_pop = SurveyPopulation(pop, survey)

    return surv_pop.frbs.fluence


def get_fluences():
    """Get fluences of populations."""
    fluences = {}

    fluences['A'] = get_local(make=MAKE, size=SIZE)
    fluences['B'] = get_further(GAMMAS[0], make=MAKE, size=SIZE)
    fluences['C'] = get_further(GAMMAS[1], make=MAKE, size=SIZE)

    return fluences


def calc_cum_hist(fluences):
    """Get the x,y for a cumulative histogram of given fluences."""
    data = {}
    for f in fluences:
        # Bin up
        fluence = fluences[f]
        min_f = np.log10(min(fluence))
        max_f = np.log10(max(fluence))
        log_bins = np.logspace(min_f, max_f, 50)
        hist, edges = np.histogram(fluence, bins=log_bins)
        n_gt_s = np.cumsum(hist[::-1])[::-1]

        x = edges[:-1]
        y = n_gt_s

        data[f] = (x, y)

    return data


def plot_logn_logs(data):
    """Plot log N log S data in a cumlative histogram."""
    fig, (ax1) = plt.subplots(1, 1)

    for key in data:
        x, y = data[key]
        ax1.step(x, y, where='post', label=key)

    plt.xlabel('Fluence (Jy ms)')
    plt.ylabel('N(>Fluence)')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim((1e-2, 1e2))
    plt.legend()
    plt.tight_layout()
    plt.savefig('plots/logn_logs_abc.pdf')


def main():
    """Run."""
    fluences = get_fluences()
    data = calc_cum_hist(fluences)
    plot_logn_logs(data)


if __name__ == '__main__':
    main()
