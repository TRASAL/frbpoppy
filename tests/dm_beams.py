"""Plot the change in DM distributions due to differing beampatterns."""
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation

from convenience import plot_aa_style, rel_path

BEAMPATTERNS = ['perfect', 'gaussian', 'airy-0', 'airy-4']
SIZE = 1e6


def get_data():
    """Get the population data."""
    # Construct population
    pop = CosmicPopulation(SIZE,
                           n_days=1,
                           name='standard_candle',
                           H_0=67.74,
                           W_m=0.3089,
                           W_v=0.6911,
                           dm_host_model='gaussian',
                           dm_host_mu=100,
                           dm_host_sigma=0,
                           dm_igm_index=1000,
                           dm_igm_sigma=None,
                           dm_mw_model='ne2001',
                           emission_range=[10e6, 10e9],
                           lum_range=[1e36, 1e36],
                           lum_index=0.,
                           n_model='sfr',
                           alpha=-1.5,
                           w_model='uniform',
                           w_range=[1., 1.],
                           w_mu=0.1,
                           w_sigma=0.,
                           si_mu=0.,
                           si_sigma=0.,
                           z_max=2.5)

    # Survey population
    pops = {}
    for b in BEAMPATTERNS:

        n_s = 0
        bp = b
        if b.startswith('airy'):
            bp, n_s = b.split('-')
            n_s = int(n_s)

        survey = Survey(name='perfect-small')
        survey.gain_pattern = bp
        survey.n_sidelobes = n_s
        surv_pop = SurveyPopulation(pop, survey)
        pops[b] = surv_pop

    return pops


def plot_dm(pops):
    """Plot resulting dispersion measure."""

    plot_aa_style()
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
    plt.savefig(rel_path(f'./plots/dm_beams.pdf'))


def main():
    pops = get_data()
    plot_dm(pops)


if __name__ == '__main__':
    main()
