"""Test the behaviour of clustered burts."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import frbpoppy.time_dists as td
from frbpoppy import pprint

from convenience import plot_aa_style, rel_path

M_BURSTS = 1100
N_DAYS = 50
N_SRCS = int(1e5)
R = 5.7
K = 0.34


def get_probabilities(normalise=False):
    """Get the number of bursts per maximum time.

    Args:
        normalise (type): Whether to normalise at each maximum time, such that
            there's a number of bursts at each time stamp having a maximum
            probability of one.

    Returns:
        array: Number of bursts per maximum time.

    """
    days = np.arange(1, N_DAYS, 1)
    m_bursts = np.arange(0, M_BURSTS, 1)
    xx, yy = np.meshgrid(m_bursts, days)
    prob = np.full((len(days), len(m_bursts)), np.nan)

    # Mask any frbs over the maximum time
    time = td.clustered(r=R, k=K, n_srcs=N_SRCS, n_days=N_DAYS, z=0)

    pprint('Masking days')
    for i, day in reversed(list(enumerate(days))):
        time[(time > day)] = np.nan
        m_bursts = (~np.isnan(time)).sum(1)
        unique, counts = np.unique(m_bursts, return_counts=True)

        try:
            prob[(i, unique)] = counts
        except IndexError:
            pprint('Please ensure M is large enough.')
            exit()

    prob = prob.T/N_SRCS

    if normalise:
        # Normalise at each maximum time (highest chance becomes one)
        prob = prob / np.nanmax(prob, axis=0)

    return prob


def plot(prob, show=False):
    """Plot the number of bursts seen over a maximum time."""
    plot_aa_style(cols=2)

    days = np.arange(1, N_DAYS, 1)
    m_bursts = np.arange(0, M_BURSTS, 1)

    prob[np.isnan(prob)] = 0

    fig, ax = plt.subplots()
    extent = [min(days)-0.5, max(days)+0.5, min(m_bursts), max(m_bursts)]
    plt.imshow(prob, norm=LogNorm(), origin='lower', aspect='auto',
               interpolation='none', extent=extent)
    plt.colorbar(orientation='vertical')
    ax.set_title('Probability of M bursts within time')
    ax.set_xlabel('Maximum time (days)')
    ax.set_ylabel('M bursts')
    plt.tight_layout()

    if show:
        plt.show()
    else:
        plt.savefig(rel_path(f'./plots/prob_m_bursts.pdf'))


if __name__ == '__main__':
    # Stops RuntimeWarnings about nan values
    np.warnings.filterwarnings('ignore')
    prob = get_probabilities(normalise=False)
    plot(prob)
