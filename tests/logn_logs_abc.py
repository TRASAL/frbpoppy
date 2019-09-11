"""Plot a log N / log S graph for three different populations."""
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import Survey, SurveyPopulation

from quick import get_cosmic_pop

MAKE = True
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

    return surv_pop.frbs.s_peak


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

    return surv_pop.frbs.s_peak


def get_s_peaks():
    """Get s_peaks of populations."""
    s_peaks = {}

    s_peaks['A'] = get_local(make=MAKE, size=SIZE)
    s_peaks['B'] = get_further(GAMMAS[0], make=MAKE, size=SIZE)
    s_peaks['C'] = get_further(GAMMAS[1], make=MAKE, size=SIZE)

    return s_peaks


def calc_cum_hist(s_peaks):
    """Get the x,y for a cumulative histogram of given s_peaks."""
    data = {}
    for s in s_peaks:
        # Bin up
        s_peak = s_peaks[s]
        min_f = np.log10(min(s_peak))
        max_f = np.log10(max(s_peak))
        log_bins = np.logspace(min_f, max_f, 50)
        hist, edges = np.histogram(s_peak, bins=log_bins)
        n_gt_s = np.cumsum(hist[::-1])[::-1]

        x = edges[:-1]
        y = n_gt_s

        data[s] = (x, y)

    return data


def plot_logn_logs(data):
    """Plot log N log S data in a cumlative histogram."""
    fig, (ax1) = plt.subplots(1, 1)

    for key in data:
        x, y = data[key]
        ax1.step(x, y, where='post', label=key)

    plt.xlabel(r'S$_{\text{min}}$ (Jy)')
    plt.ylabel(r'N(>S$_{\text{min}}$)')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim((1e-3, 1e1))
    plt.legend()
    plt.tight_layout()
    plt.savefig('plots/logn_logs_abc.pdf')


def main():
    """Run."""
    s_peaks = get_s_peaks()
    data = calc_cum_hist(s_peaks)
    plot_logn_logs(data)


if __name__ == '__main__':
    main()
