"""Use simple rate comparisions, try predicting event rates."""
import numpy as np
import matplotlib.pyplot as plt
import os

from frbpoppy import Survey

from convenience import plot_aa_style, rel_path

ALPHAS = np.around(np.linspace(-0.2, -2.5, 7), decimals=2)
SURVEYS = ('palfa', 'htru', 'askap-fly')


def compare_surveys(surv1, surv2, alpha):
    """Event rate surv1 / Event rate surv2 for an alpha."""
    omega = surv1.beam_size_fwhm/surv2.beam_size_fwhm
    T_rec = surv1.T_rec/surv2.T_rec
    gain = surv1.gain/surv2.gain
    beta = surv1.beta/surv2.beta
    SEFD = T_rec*beta/gain
    bw = surv1.bw/surv2.bw
    S_min = surv1.snr_limit/surv2.snr_limit

    return omega * (SEFD * S_min)**alpha * (bw)**(-alpha/2)


def toy_rates(surveys=SURVEYS, alphas=ALPHAS):
    """Use a toy model to scale detection rates to various alphas."""
    rates = {}

    for surv in surveys:

        # Get survey parameters
        surv1 = Survey(surv, beam_pattern='perfect', n_sidelobes=0.5)
        surv2 = Survey('htru', beam_pattern='perfect', n_sidelobes=0.5)

        # Calculate rate per alpha
        rate = []
        for alpha in alphas:
            rate.append(compare_surveys(surv1, surv2, alpha))
        rates[surv] = rate

    return rates


def main():
    """Plot toy rates."""
    rates = toy_rates()

    for surv in rates:
        rate = rates[surv]

        plot_aa_style()

        plt.plot(ALPHAS, rate, label=surv)

        plt.xlabel(r'$\alpha_{\text{in}}$')
        plt.ylabel(r'Events / htru')
        plt.xlim((min(ALPHAS), max(ALPHAS)))
        plt.yscale('log')
        plt.legend()
        plt.grid()
        plt.gca().invert_xaxis()
        plt.tight_layout()
        plt.savefig(rel_path('./plots/toy_rates.pdf'))


if __name__ == '__main__':
    main()
