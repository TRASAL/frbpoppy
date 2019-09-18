# -*- coding: future_fstrings -*-
"""Plot the change in DM distributions due to differing beampatterns."""
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import Survey, SurveyPopulation

from quick import get_cosmic_pop

MAKE = False
BEAMPATTERNS = ['perfect', 'gaussian', 'airy-0', 'airy-4']
SIZE = 'medium'


def do_survey():

    # Construct population
    pop = get_cosmic_pop('standard_candle',
                         SIZE,
                         load=True,
                         overwrite=MAKE)

    # Survey population
    pops = {}
    for b in BEAMPATTERNS:

        n_s = 0
        bp = b
        if b.startswith('airy'):
            bp, n_s = b.split('-')
            n_s = int(n_s)

        survey = Survey(name='perfect-small', gain_pattern=bp, n_sidelobes=n_s)
        surv_pop = SurveyPopulation(pop, survey)
        pops[b] = surv_pop

    return pops


def plot_dm(pops):

    f, (ax1) = plt.subplots(1, 1)

    for i, beam in enumerate(pops):
        pop = pops[beam]
        dm = pop.frbs.dm
        s_peak = pop.frbs.s_peak

        limit = 1e-10
        dm = dm[(s_peak > limit)]

        print(f'{len(dm)} FRBs in graph of {pop.name}')

        n, bin_edges = np.histogram(dm, bins=75)
        n = n/max(n)
        bins = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Add edges to histogram
        bin_diff = np.diff(bins)[-1]
        bins = np.insert(bins, 0, bins[0] - bin_diff)
        bins = np.insert(bins, len(bins), bins[-1] + bin_diff)
        n = np.insert(n, 0, 0)
        n = np.insert(n, len(n), 0)

        ax1.step(bins, n, where='mid', label=beam)

    ax1.set_xlabel(r'DM (pc cm$^{-3}$)')
    ax1.set_xlim([0, 3000])
    ax1.set_ylabel(r'Fraction')
    ax1.set_ylim([0, 1])
    ax1.legend()

    plt.tight_layout()
    plt.savefig(f'plots/dm_beams.pdf')


def main():
    pops = do_survey()
    plot_dm(pops)


if __name__ == '__main__':
    main()
