"""Plot corner plots showing best fit for rates.

Linked with the ideas in cube.py -> generating a range of parameters.

TODO: In progress
"""
import numpy as np
import os
from glob import glob
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from frbpoppy import paths, unpickle, poisson_interval

from tests.convenience import plot_aa_style, rel_path
from alpha_real import EXPECTED

# # Set parameters needed for generating a rate cube
# GENERATE = False
# PLOT = True
vs = {'alpha': np.linspace(-1, -2, 5)[::-1][1:4],
      'li': np.linspace(-2, 0, 5),  # Luminosity function index
      'si': np.linspace(-2, 2, 5)}  # Spectral index
SURVEYS = ('askap-fly', 'fast', 'htru', 'apertif', 'palfa')


def get_pops(alpha='*', li='*', si='*', survey='*'):
    filename = f'complex_alpha_{alpha}_lum_{li}_si_{si}_{survey}.p'
    filter = os.path.join(paths.populations(), filename)
    pop_paths = glob(filter)

    pops = []
    for path in pop_paths:
        if '_for_plotting' not in path:
            pops.append(unpickle(path))
    return pops


def make_mesh(x_par, y_par, survey):
    v = np.zeros([len(vs[x_par]), len(vs[y_par])])
    print(x_par, y_par)
    for i, x_val in enumerate(vs[x_par]):
        for j, y_val in enumerate(vs[y_par]):
            print(x_val, y_val)

            pops = get_pops(survey=survey, **{x_par: x_val, y_par: y_val})

            if pops:
                rates = []
                for pop in pops:
                    rate = pop.source_rate.det / pop.source_rate.days
                    rate_err = poisson_interval(pop.source_rate.det, sigma=1)
                    if pop.source_rate.det == 0:
                        rate = np.nan
                    rates.append(rate)

                mean_rate = np.nanmean(rates)
                exp_rate = EXPECTED[survey][0]/EXPECTED[survey][1]
                v[i, j] = np.abs(mean_rate - exp_rate)
    return v.T


for survey in SURVEYS:

    plot_aa_style()
    plt.rcParams["figure.figsize"] = (5.75373, 5.75373)

    fig, axes = plt.subplots(3, 3)

    args = {'cmap': 'viridis', 'norm': LogNorm(), 'vmin': 1e-4, 'vmax': 1e6}
    axes[2, 0].imshow(make_mesh('alpha', 'si', survey), **args)
    axes[1, 0].imshow(make_mesh('alpha', 'li', survey), **args)
    im = axes[2, 1].imshow(make_mesh('li', 'si', survey), **args)
    fig.colorbar(im, ax=axes[0, 2])
    axes[1, 1].set_title(survey)

    axes[2, 0].set_xlabel('alpha')
    axes[2, 0].set_ylabel('si')
    axes[2, 1].set_xlabel('li')
    axes[2, 2].set_xlabel('si')
    axes[1, 0].set_ylabel('li')
    axes[0, 0].set_ylabel('alpha')
    axes[0, 1].set_axis_off()
    axes[0, 2].set_axis_off()
    axes[1, 2].set_axis_off()
    axes[0, 0].set_yticks([])
    axes[1, 1].set_yticks([])
    axes[2, 2].set_yticks([])

    # Set up axes
    for ax in axes.flat:
        ax.label_outer()
        x_par = ax.get_xlabel()
        if x_par:
            ax.set_xticks(np.arange(len(vs[x_par])))
            ax.set_xticklabels(vs[x_par])
        y_par = ax.get_ylabel()
        if y_par:
            ax.set_yticks(np.arange(len(vs[y_par])))
            ax.set_yticklabels(vs[y_par])

    plt.tight_layout()
    plt.savefig(rel_path(f'./plots/corner_{survey}.pdf'))
