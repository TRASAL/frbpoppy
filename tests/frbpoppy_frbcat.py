"""Plot the DM distribution obtained with frbpoppy against frbcat results."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from scipy.stats import ks_2samp

from frbpoppy import Survey, Frbcat, pprint

from quick import get_cosmic_pop, get_survey_pop

MAKE = False
OBSERVE = False
PLOT = True
SIZE = 'large'
TELESCOPES = ['parkes', 'askap']


def hist(parameter, bins='lin', n_bins=25):
    """Bin up a parameter either in a lin or log space.

    Why is this not a standard option in numpy or matplotlib?

    Args:
        parameter (array): To be binned
        bins (str): Either 'lin' or 'log'
        n_bins (int): Number of bins. Can be overriden internally

    Returns:
        tuple: bin centers, values per bin

    """
    # Drop NaN-values
    parameter = parameter[~np.isnan(parameter)]

    # Determine number of bins
    if len(parameter) < 50:
        n_bins = 15
    if len(parameter) > 500:
        n_bins = 50

    # Determine type of binning
    if bins == 'lin':
        bins = n_bins
    elif bins == 'log':
        min_f = np.log10(min(parameter))
        max_f = np.log10(max(parameter))
        bins = np.logspace(min_f, max_f, n_bins)

    # Bin
    n, bins = np.histogram(parameter, bins=bins)
    n = n/max(n)  # Normalise
    bin_centres = (bins[:-1] + bins[1:]) / 2

    return bin_centres, n


def plot_dists(surv_pop, telescope):
    """Plot the fluence and DM distribution of a surveyed population.

    Args:
        surv_pop (Population): Population from which to plot
        telescope (str): Name of the telescope with which to compare the
            distribution. Necessary for Frbcat.

    """
    # Use a nice font for axes
    plt.rc('text', usetex=True)

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
    ax2.step(*hist(fluence_frbpoppy, bins='log'), where='mid',
             label='frbpoppy', linestyle='dashed')

    fluence_frbcat = df[df.telescope == telescope].fluence
    ax2.step(*hist(fluence_frbcat, bins='log'), where='mid', label='frbcat')

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


def main():
    """Run main part of code."""
    # Be a bit smart about which populations need to be loaded
    load = True
    if not MAKE and not OBSERVE:
        load = False

    pop = get_cosmic_pop('standard', SIZE, load=load, overwrite=MAKE)

    for telescope in TELESCOPES:

        pattern = 'airy'
        if telescope == 'parkes':
            pattern = telescope

        s = telescope
        if telescope == 'askap':
            s = 'askap-fly'

        survey = Survey(s, gain_pattern=pattern, n_sidelobes=1)

        surv_pop = get_survey_pop(pop, survey, overwrite=OBSERVE)

        if PLOT:
            plot_dists(surv_pop, telescope)


if __name__ == '__main__':
    main()
