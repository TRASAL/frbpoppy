"""Compare rate calculations per alpha for the two askap settings."""
import numpy as np
import matplotlib.pyplot as plt

from tests.convenience import plot_aa_style, rel_path
from alpha_complex import complex_rates

REMAKE = False
SIZE = 1e4
SURVEYS = ['htru', 'askap-fly', 'askap-incoh']
ALPHAS = np.around(np.linspace(-0.2, -2.5, 7), decimals=2)


def main():
    """Get detection rates for surveys."""
    complex = complex_rates(remake=REMAKE,
                            alphas=ALPHAS,
                            size=SIZE,
                            surveys=SURVEYS)

    # Plot population event rates
    plot_rates(complex)


def plot_rates(rates):
    """Plot detection rates for askap surveys."""
    plot_aa_style()

    fig, (ax1) = plt.subplots(1, 1)
    cmap = plt.get_cmap('tab10')
    ax1.set_xlim((min(ALPHAS)+.1, max(ALPHAS)-.1))
    ax1.set_yscale('log', nonposy='mask')

    # Plot complex versus toy
    for i, surv in enumerate(SURVEYS):
        ax1.plot(ALPHAS, rates[surv], color=cmap(i+1), label=surv,
                 linestyle='dashed')

    # Plot layout options
    # Set up axes
    ax1.set_xlabel(r'$\alpha_{\text{in}}$')
    ax1.invert_xaxis()
    ax1.set_ylabel('Events / htru')

    plt.legend()
    plt.savefig(rel_path('plots/askap_rates.pdf'), bbox_inches='tight')


if __name__ == '__main__':
    main()
