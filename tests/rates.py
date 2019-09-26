"""Show frbpoppy matches analytical models and predict the event rates."""
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import os

from rates_toy import toy_rates
from rates_real import real_rates
from rates_simple import simple_rates
from rates_complex import complex_rates

REMAKE = True
SIZE = 1e8
SURVEYS = ('palfa', 'htru', 'askap-fly', 'askap-incoh')
ALPHAS = np.around(np.linspace(-0.2, -2.5, 7), decimals=2)


def plot(toy, simple, complex, real):
    """Plot rates panel."""

    surveys = SURVEYS[:-1]

    # Use a nice font for axes
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    plt.style.use('./aa.mplstyle')  # Use A&A styling for plots
    plt.rcParams["figure.figsize"] = (5.75373, 3.556)

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    cmap = plt.get_cmap('tab10')
    ax1.set_xlim((min(ALPHAS)+.1, max(ALPHAS)-.1))
    ax2.set_xlim((min(ALPHAS)+.1, max(ALPHAS)-.1))
    ax1.set_yscale('log', nonposy='mask')
    ax2.set_yscale('log', nonposy='mask')

    # Plot simple versus toy
    for i, surv in enumerate(surveys):
        ax1.plot(ALPHAS, toy[surv], color=cmap(i), linestyle='dotted',
                 zorder=0)
        ax1.plot(ALPHAS, simple[surv], zorder=1)

    # Plot complex expectations
    for i, surv in enumerate(surveys):
        ax2.plot(ALPHAS, complex[surv], color=cmap(i), linestyle='dashed',
                 zorder=1)

    # Plot real event rate boxes
    ma, mi = ax2.get_xlim()
    ma -= 0.05
    mi += 0.05
    size = 0.13
    z = 0
    for i, surv in enumerate(surveys):

        central, min_r, max_r = real[surv]

        left = mi - size
        right = ma + size

        x, y = zip(*[(ma, max_r), (right, max_r), (right, min_r), (ma, min_r)])
        ax1.fill(x, y, color=cmap(i), zorder=z)
        ax1.plot([ma, right+0.08], [central, central], color=cmap(i), zorder=z)

        x, y = zip(*[(mi, max_r), (left, max_r), (left, min_r), (mi, min_r)])
        ax2.fill(x, y, color=cmap(i), zorder=z)
        ax2.plot([mi, left-0.08], [central, central], color=cmap(i), zorder=z)

        size -= 0.02
        z += 1

    # Plot layout options
    # Set up axes
    ax1.set_xlabel(r'$\alpha_{\text{in}}$')
    ax1.invert_xaxis()
    ax1.set_ylabel('Events / htru')
    ax1.yaxis.set_ticks_position('left')
    ax1.title.set_text(r'\textit{Simple} populations')

    ax2.set_xlabel(r'$\alpha_{\text{in}}$')
    ax2.invert_xaxis()
    ax2.yaxis.set_ticks_position('right')
    ax2.tick_params(labelright=False)
    ax2.title.set_text(r'\textit{Complex} populations')

    # Set up layout options
    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)

    # Add legend elements
    elements = []
    for i, surv in enumerate(surveys):
        c = cmap(i)
        line = Line2D([0], [0], color=c)
        label = surv
        elements.append((line, label))

    # Add gap in legend
    elements.append((Line2D([0], [0], color='white'), ''))

    # Add line styles
    n = 'analytical'
    elements.append((Line2D([0], [0], color='gray', linestyle='dotted'), n))
    elements.append((Line2D([0], [0], color='gray'), 'simple'))
    elements.append((Line2D([0], [0], color='gray', linestyle='dashed'),
                     'complex'))

    # Add gap in legend
    elements.append((Line2D([0], [0], color='white'), ''))

    elements.append((Patch(facecolor='gray', edgecolor='gray', alpha=0.6),
                     'real'))

    lines, labels = zip(*elements)
    plt.legend(lines, labels, bbox_to_anchor=(1.04, 0.5), loc="center left")

    plt.savefig('plots/rates.pdf', bbox_inches='tight')


def main():
    toy = toy_rates(surveys=SURVEYS,
                    alphas=ALPHAS)

    simple = simple_rates(remake=REMAKE,
                          alphas=ALPHAS,
                          size=SIZE,
                          surveys=SURVEYS)

    complex = complex_rates(remake=REMAKE,
                            alphas=ALPHAS,
                            size=SIZE,
                            surveys=SURVEYS)

    real = real_rates(surveys=SURVEYS)

    plot(toy, simple, complex, real)


if __name__ == '__main__':
    main()
