"""Show frbpoppy matches analytical models and predict the event rates."""
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

from rates_toy import toy_rates
from rates_real import real_rates
from rates_simple import simple_rates
from rates_complex import complex_rates

MAKE = False
OBSERVE = False
SIZE = 'large'
SURVEYS = ('palfa', 'htru', 'askap-fly')
ALPHAS = np.around(np.linspace(-0.2, -2.5, 7), decimals=2)

james = [-2.2-0.47, -2.2+0.47]  # James, Ekers et al. 2018


def plot(toy, simple, complex, real):

    fig, (ax1) = plt.subplots(1, 1, sharey=True)
    cmap = plt.get_cmap('tab10')
    ax1.set_xlim((min(ALPHAS)+.1, max(ALPHAS)-.1))
    ax1.set_yscale('log', nonposy='mask')

    # Plot simple versus toy
    # for i, surv in enumerate(SURVEYS):
    #     ax1.plot(ALPHAS, toy[surv], color=cmap(i), linestyle='dotted',
    #              zorder=0)
    #     ax1.plot(ALPHAS, simple[surv], zorder=1)

    # Plot complex expectations
    for i, surv in enumerate(SURVEYS):
        ax1.plot(ALPHAS, complex[surv], color=cmap(i), linestyle='dashed',
                 zorder=1)

    # Plot real event rate boxes
    ma, mi = ax1.get_xlim()
    ma -= 0.05
    mi += 0.05
    size = 0.13
    z = 0
    for i, surv in enumerate(SURVEYS):

        central, min_r, max_r = real[surv]

        left = mi - size
        right = ma + size

        x, y = zip(*[(ma, max_r), (right, max_r), (right, min_r), (ma, min_r)])
        ax1.fill(x, y, color=cmap(i), zorder=z)
        ax1.plot([ma, right+0.08], [central, central], color=cmap(i), zorder=z)

        size -= 0.02
        z += 1

    # Plot layout options
    # Set up axes
    ax1.set_xlabel(r'$\alpha$')
    ax1.invert_xaxis()
    ax1.set_ylabel('Events / htru')
    ax1.yaxis.set_ticks_position('left')

    # Add legend elements
    elements = []
    for i, surv in enumerate(SURVEYS):
        c = cmap(i)
        line = Line2D([0], [0], color=c)
        label = surv
        elements.append((line, label))

    # Add gap in legend
    elements.append((Line2D([0], [0], color='white'), ''))

    # Add line styles
    n = 'analytical'
    # elements.append((Line2D([0], [0], color='gray', linestyle='dotted'), n))
    elements.append((Line2D([0], [0], color='gray', linestyle='dashed'),
                     'complex'))
    elements.append((Patch(facecolor='gray', edgecolor='gray', alpha=0.6),
                     'real'))
    # elements.append((Line2D([0], [0], color='gray'), 'simple'))
    # elements.append((Patch(facecolor='white', edgecolor='gray', hatch='+',
    #                        linewidth=0.1), 'prediction'))

    lines, labels = zip(*elements)
    plt.legend(lines, labels, bbox_to_anchor=(1.04, 0.5), loc="center left")

    plt.savefig('plots/rates_presentation.pdf', bbox_inches='tight')


def main():
    toy = toy_rates(surveys=SURVEYS,
                    alphas=ALPHAS)

    simple = simple_rates(make=MAKE,
                          observe=OBSERVE,
                          alphas=ALPHAS,
                          size=SIZE,
                          surveys=SURVEYS)

    complex = complex_rates(make=MAKE,
                            observe=OBSERVE,
                            alphas=ALPHAS,
                            size=SIZE,
                            surveys=SURVEYS)

    real = real_rates(surveys=SURVEYS)

    plot(toy, simple, complex, real)


if __name__ == '__main__':
    main()
