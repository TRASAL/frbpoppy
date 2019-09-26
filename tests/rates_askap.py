"""Compare rate calculations per alpha for the two askap settings."""
import numpy as np
import matplotlib.pyplot as plt
import os

from rates_complex import complex_rates

REMAKE = False
SIZE = 1e8
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
    # Change working directory
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    # Use A&A styling for plotting
    plt.style.use('./aa.mplstyle')

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
    ax1.set_xlabel(r'$\alpha$')
    ax1.invert_xaxis()
    ax1.set_ylabel('Events / htru')

    plt.legend()
    plt.savefig('plots/askap_rates.pdf', bbox_inches='tight')


if __name__ == '__main__':
    main()
