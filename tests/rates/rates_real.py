"""Calculate the real frb detection rates."""
from scipy.stats import chi2, norm
from scipy.integrate import quad
import matplotlib.pyplot as plt
import numpy as np

from tests.convenience import plot_aa_style, rel_path

EXPECTED = {'htru': [9, 0.551 / 1549 / 24],  # N_frbs, N_days
            'apertif': [9, 1100/24],  # 1100 hours
            'askap-fly': [20, 32840 * 8 / 24],
            'palfa': [1, 24.1],
            'guppi': [0.4, 81],  # 0.4 is my own assumption
            'fast': [1, 1500/24]
            }

SURVEYS = ('askap-fly', 'fast')
ALPHAS = np.around(np.linspace(-0.5, -2.0, 7), decimals=2)


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


def real_rates(surveys=SURVEYS):
    """Calculate the EXPECTED rates (all scaled to a survey)."""
    rates = {}
    scale_to = surveys[0]

    for surv in surveys:

        if surv not in EXPECTED:
            continue

        # Plot EXPECTED rate
        exp_n = EXPECTED[surv][0]
        exp_days = EXPECTED[surv][1]
        scale_n = EXPECTED[scale_to][0]
        scale_days = EXPECTED[scale_to][1]

        exp_min, exp_max = poisson_interval(exp_n, sigma=2)

        exp = (exp_n/exp_days) / (scale_n/scale_days)
        exp_min *= (1/exp_days) / (scale_n/scale_days)
        exp_max *= (1/exp_days) / (scale_n/scale_days)

        rates[surv] = (exp, exp_min, exp_max)

    return rates


def main():
    """Plot real rate regions per alpha."""
    plot_aa_style()

    rates = real_rates()

    for surv in rates:
        middle, top, bottom = rates[surv]
        left = min(ALPHAS)
        right = max(ALPHAS)
        x = [left, right, right, left]
        y = [top, top, bottom, bottom]
        plt.fill(x, y, alpha=0.25)
        plt.plot([left, right], [middle]*2, label=surv, linestyle='dashed')

    plt.xlabel(r'$\alpha_{\text{in}}$')
    plt.ylabel(rf'Events / {SURVEYS[0]}')
    plt.yscale('log')
    plt.legend()
    plt.gca().invert_xaxis()
    plt.tight_layout()
    plt.savefig(rel_path('./plots/real_rates.pdf'))


if __name__ == '__main__':
    main()
