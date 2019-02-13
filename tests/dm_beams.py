"""Plot the change in DM distributions due to differing beampatterns."""
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import Survey, SurveyPopulation

from quick import get_cosmic_pop

MAKE = False
BEAMPATTERNS = ['perfect', 'gaussian', 'airy-0', 'airy-8']
SIZE = 'small'


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

    for beam in pops:
        pop = pops[beam]
        dm = pop.frbs.dm
        s_peak = pop.frbs.s_peak

        limit = 1e-9
        dm = dm[(s_peak > limit)]

        print(f'{len(dm)} FRBs in graph of {pop.name}')

        n, bins = np.histogram(dm, bins=50)
        n = n/max(n)
        bincentres = (bins[:-1] + bins[1:]) / 2
        ax1.step(bincentres, n, where='mid', label=beam)

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
