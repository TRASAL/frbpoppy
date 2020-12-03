"""Show frbpoppy matches analytical models and predict the event rates."""
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np

from tests.convenience import plot_aa_style, rel_path
from alpha_complex import complex_rates
from alpha_real import real_rates
from alpha_simple import simple_rates
from alpha_analytical import analytical_rates

REMAKE = True
SIZE = 1e4
SURVEYS = ('askap-fly', 'fast-crafts', 'parkes-htru', 'wsrt-apertif', 'arecibo-palfa')
ELEMENTS = {'analytical': True, 'real': True, 'simple': False, 'complex': True}
ALPHAS = np.around(np.linspace(-0.5, -2.0, 7), decimals=2)


def plot(analytical=True, simple=False, complex=False, real=True):
    """Plot rates panel."""
    surveys = SURVEYS

    # If needing two panels
    if simple and complex:
        plot_aa_style(cols=2)
        fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    else:
        plot_aa_style(cols=1)
        fig, (ax1) = plt.subplots(1, 1)
        ax2 = ax1

    cmap = plt.get_cmap('tab10')
    ax1.set_xlim((min(ALPHAS)+.1, max(ALPHAS)-.1))
    ax2.set_xlim((min(ALPHAS)+.1, max(ALPHAS)-.1))
    ax1.set_yscale('log', nonposy='mask')
    ax2.set_yscale('log', nonposy='mask')

    # Plot simple versus analytical
    for i, surv in enumerate(surveys):
        if analytical:
            ax1.plot(ALPHAS, analytical[surv], color=cmap(i),
                     linestyle='dotted', zorder=0)
        if simple:
            ax1.plot(ALPHAS, simple[surv], zorder=1)

    # Plot complex expectations
    for i, surv in enumerate(surveys):
        if complex:
            ax2.plot(ALPHAS, complex[surv], color=cmap(i), linestyle='dashed',
                     zorder=1)

    # Plot real event rate boxes
    if real:
        ma, mi = ax2.get_xlim()
        ma -= 0.05
        mi += 0.05
        size = len(surveys)*0.015+0.05  # 0.13
        z = 0
        for i, surv in enumerate(surveys):

            central, min_r, max_r = real[surv]

            left = mi - size
            right = ma + size

            x, y = zip(*[(ma, max_r), (right, max_r), (right, min_r),
                         (ma, min_r)])
            ax1.fill(x, y, color=cmap(i), zorder=z)
            ax1.plot([ma, right+0.08], [central, central], color=cmap(i),
                     zorder=z)

            x, y = zip(*[(mi, max_r), (left, max_r), (left, min_r),
                         (mi, min_r)])
            ax2.fill(x, y, color=cmap(i), zorder=z)
            ax2.plot([mi, left-0.08], [central, central], color=cmap(i),
                     zorder=z)

            size -= 0.015
            z += 1

    # Plot layout options
    # Set up axes
    ax1.set_xlabel(r'$\alpha_{\text{in}}$')
    ax1.invert_xaxis()
    ax1.set_ylabel(f'Events / {SURVEYS[0]}')
    ax1.yaxis.set_ticks_position('left')
    if simple:
        ax1.title.set_text(r'\textit{Simple} populations')

    if simple and complex:
        ax2.set_xlabel(r'$\alpha_{\text{in}}$')
        ax2.invert_xaxis()
        ax2.yaxis.set_ticks_position('right')
        ax2.tick_params(labelright=False)
        if complex:
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
    if analytical:
        n = 'analytical'
        elements.append((Line2D([0], [0], color='gray', linestyle='dotted'),
                        n))
    if simple:
        elements.append((Line2D([0], [0], color='gray'), 'simple'))
    elements.append((Line2D([0], [0], color='gray', linestyle='dashed'),
                     'complex'))

    # Add gap in legend
    elements.append((Line2D([0], [0], color='white'), ''))

    if real:
        elements.append((Patch(facecolor='gray', edgecolor='gray', alpha=0.6),
                         'real'))

    lines, labels = zip(*elements)
    plt.legend(lines, labels, bbox_to_anchor=(1.04, 0.5), loc="center left")

    plt.savefig(rel_path('plots/rates_overview.pdf'), bbox_inches='tight')


def main():
    """Get rates."""

    if ELEMENTS['analytical']:
        ELEMENTS['analytical'] = analytical_rates(surveys=SURVEYS,
                                                  alphas=ALPHAS)

    if ELEMENTS['simple']:
        ELEMENTS['simple'] = simple_rates(remake=REMAKE,
                                          alphas=ALPHAS,
                                          size=SIZE,
                                          surveys=SURVEYS)

    if ELEMENTS['complex']:
        ELEMENTS['complex'] = complex_rates(remake=REMAKE,
                                            alphas=ALPHAS,
                                            size=SIZE,
                                            surveys=SURVEYS)

    if ELEMENTS['real']:
        ELEMENTS['real'] = real_rates(surveys=SURVEYS)

    plot(**ELEMENTS)


if __name__ == '__main__':
    main()
