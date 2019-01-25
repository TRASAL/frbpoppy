"""Plot the DM distribution obtained with frbpoppy against frbcat results."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from scipy.stats import ks_2samp

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, unpickle
from frbpoppy import Frbcat, pprint

from gen_standard import gen_standard

MAKE = False
OBSERVE = False


def hist(parameter, bins='lin', n_bins=25):

    # Drop NaN-values
    parameter = parameter[~np.isnan(parameter)]

    n_bins = 30
    if len(parameter) < 50:
        n_bins = 15

    if bins == 'lin':
        bins = n_bins
    elif bins == 'log':
        min_f = np.log10(min(parameter))
        max_f = np.log10(max(parameter))
        bins = np.logspace(min_f, max_f, n_bins)

    n, bins = np.histogram(parameter, bins=bins)
    n = n/max(n)
    bin_centres = (bins[:-1] + bins[1:]) / 2
    return bin_centres, n

def plot_dists(surv_pop, telescope):
    # Plot dm distribution
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

    dm_frbpoppy = surv_pop.frbs.dm
    pprint(f'Number of detected FRBs: {len(dm_frbpoppy)}')
    bin_centres, n = hist(dm_frbpoppy)
    ax1.step(bin_centres, n, where='mid')

    df = Frbcat().df
    dm_frbcat = df[df.telescope==telescope].dm
    bin_centres, n = hist(dm_frbcat)
    ax1.step(bin_centres, n, where='mid')

    # Compare distributions
    ks = ks_2samp(dm_frbpoppy, dm_frbcat)
    text = fr'$p={round(ks[1], 2)}$'
    anchored_text = AnchoredText(text, loc=1, borderpad=1., frameon=False)
    ax1.add_artist(anchored_text)

    ax1.set_xlabel(r'DM ($\textrm{pc}\ \textrm{cm}^{-3}$)')
    ax1.set_ylabel('Fraction')
    ax1.set_ylim([0, 1.1])
    ax1.set_xlim([0, 3500])

    plt.tight_layout()

    # Plot fluence distributions

    fluence_frbpoppy = surv_pop.frbs.fluence
    bin_centres, n = hist(fluence_frbpoppy, bins='log')
    ax2.step(bin_centres, n, where='mid', label='frbpoppy')

    fluence_frbcat = df[df.telescope==telescope].fluence
    bin_centres, n = hist(fluence_frbcat, bins='log')
    ax2.step(bin_centres, n, where='mid', label='frbcat')

    # Compare distributions
    ks = ks_2samp(fluence_frbpoppy, fluence_frbcat)
    text = fr'$p={round(ks[1], 2)}$'
    anchored_text = AnchoredText(text, loc=1, borderpad=1., frameon=False)
    ax2.add_artist(anchored_text)

    ax2.set_xlabel(r'Fluence (Jy ms)')
    ax2.set_ylim([0, 1.1])
    ax2.set_xlim([5e-1,1e4])
    plt.xscale('log')

    plt.figlegend(loc='upper center', ncol=2, framealpha=1)

    plt.tight_layout()
    plt.savefig(f'plots/frbpoppy_{telescope}.pdf')

def main():

    if MAKE:
        pop = gen_standard()

    if OBSERVE:

        if not MAKE:
            pop = unpickle('standard')

        survey = Survey('parkes', gain_pattern='parkes')
        surv_pop = SurveyPopulation(pop, survey)
        surv_pop.name = 'standard_parkes'
        surv_pop.save()

    else:
        surv_pop = unpickle('standard_parkes')

    plot_dists(surv_pop, 'parkes')

if __name__ == '__main__':
    main()
