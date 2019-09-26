"""Plot the DM distribution obtained with frbpoppy against frbcat results."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from scipy.stats import ks_2samp

from frbpoppy import CosmicPopulation, Survey, Frbcat, pprint, LargePopulation
from frbpoppy import paths, unpickle

from easy_hist import hist

PLOT = True
SIZE = 1e8
REMAKE = True
TELESCOPES = ['parkes', 'askap']


def plot_dists(surv_pop, telescope):
    """Plot the fluence and DM distribution of a surveyed population.

    Args:
        surv_pop (Population): Population from which to plot
        telescope (str): Name of the telescope with which to compare the
            distribution. Necessary for Frbcat.

    """
    # Use a nice font for axes
    plt.rc('text', usetex=True)
    plt.rcParams["figure.figsize"] = (5.75373, 3.556)

    # Plot dm distribution
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

    dm_frbpoppy = surv_pop.frbs.dm
    pprint(f'Number of detected FRBs: {len(dm_frbpoppy)}')
    ax1.step(*hist(dm_frbpoppy), where='mid', linestyle='dashed')

    df = Frbcat().df
    dm_frbcat = df[df.telescope == telescope].dm
    ax1.step(*hist(dm_frbcat), where='mid')

    # Compare distributions
    ks = ks_2samp(dm_frbpoppy, dm_frbcat)
    text = fr'$p={round(ks[1], 2)}$'
    anchored_text = AnchoredText(text, loc=1, borderpad=1., frameon=False)
    ax1.add_artist(anchored_text)

    ax1.set_xlabel(r'DM ($\textrm{pc}\ \textrm{cm}^{-3}$)')
    ax1.set_ylabel('Fraction')
    ax1.set_ylim([0, 1.1])
    ax1.set_xlim([0, 3500])

    # Plot fluence distributions
    fluence_frbpoppy = surv_pop.frbs.fluence
    ax2.step(*hist(fluence_frbpoppy, bin_type='log'), where='mid',
             label='frbpoppy', linestyle='dashed')

    fluence_frbcat = df[df.telescope == telescope].fluence
    ax2.step(*hist(fluence_frbcat, bin_type='log'), where='mid',
             label='frbcat')

    # Compare distributions
    ks = ks_2samp(fluence_frbpoppy, fluence_frbcat)
    text = fr'$p={round(ks[1], 2)}$'
    anchored_text = AnchoredText(text, loc=1, borderpad=1., frameon=False)
    ax2.add_artist(anchored_text)

    ax2.set_xlabel(r'Fluence (Jy ms)')
    ax2.set_ylim([0, 1.1])
    ax2.set_xlim([5e-1, 1e4])
    plt.xscale('log')

    plt.figlegend(loc='upper center', ncol=2, framealpha=1)

    plt.tight_layout()
    plt.savefig(f'plots/frbpoppy_{telescope}.pdf')
    plt.clf()


def get_data():
    """Get survey populations."""
    # Don't always regenerate a population
    if REMAKE is False:
        # Check where a possible population would be located
        path = ''
        surv_pops = []
        for telescope in TELESCOPES:
            if telescope == 'askap':
                telescope = 'askap-fly'
            name = f'{telescope}'
            path = paths.populations() + name + '.p'
            surv_pops.append(unpickle(path))

        return surv_pops

    cosmic_pop = CosmicPopulation.complex(SIZE, generate=False)

    surveys = []
    for telescope in TELESCOPES:

        pattern = 'airy'
        if telescope == 'parkes':
            pattern = telescope

        s = telescope
        if telescope == 'askap':
            s = 'askap-fly'

        surveys.append(Survey(s, gain_pattern=pattern, n_sidelobes=1))

    return LargePopulation(cosmic_pop, *surveys).pops


def main():
    """Run main part of code."""

    surv_pops = get_data()

    if PLOT:
        for i, telescope in enumerate(TELESCOPES):
            plot_dists(surv_pops[i], telescope)


if __name__ == '__main__':
    main()
