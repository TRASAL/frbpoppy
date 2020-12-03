"""Use simple rate comparisions, try predicting event rates."""
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import Survey

from tests.convenience import plot_aa_style, rel_path

ALPHAS = np.around(np.linspace(-0.2, -2.5, 7), decimals=2)
SURVEYS = ('parkes-htru', 'arecibo-palfa', 'askap-fly', 'fast-crafts')


def compare_surveys(surv1, surv2, alpha):
    """Event rate surv1 / Event rate surv2 for an alpha."""
    omega = surv1.beam_size_at_fwhm/surv2.beam_size_at_fwhm
    T_rec = surv1.T_rec/surv2.T_rec
    gain = surv1.gain/surv2.gain
    beta = surv1.beta/surv2.beta
    SEFD = T_rec*beta/gain
    bw = surv1.bw/surv2.bw
    S_min = surv1.snr_limit/surv2.snr_limit

    return omega * (SEFD * S_min)**alpha * (bw)**(-alpha/2)


def analytical_rates(surveys=SURVEYS, alphas=ALPHAS):
    """Use a analytical model to scale detection rates to various alphas."""
    rates = {}

    for surv in surveys:

        # Get survey parameters
        surv1 = Survey(surv)
        surv1.set_beam('perfect', n_sidelobes=0.5)
        surv2 = Survey(surveys[0])
        surv2.set_beam('perfect', n_sidelobes=0.5)

        # Calculate rate per alpha
        rate = []
        for alpha in alphas:
            rate.append(compare_surveys(surv1, surv2, alpha))
        rates[surv] = rate

    return rates


def main():
    """Plot analytical rates."""
    rates = analytical_rates()

    for surv in rates:
        rate = rates[surv]

        plot_aa_style()

        plt.plot(ALPHAS, rate, label=surv)

        plt.xlabel(r'$\alpha_{\text{in}}$')
        plt.ylabel(rf'Events / {SURVEYS[0]}')
        plt.xlim((min(ALPHAS), max(ALPHAS)))
        plt.yscale('log')
        plt.legend()
        plt.gca().invert_xaxis()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(rel_path('./plots/rates_analytical.pdf'))


if __name__ == '__main__':
    main()
