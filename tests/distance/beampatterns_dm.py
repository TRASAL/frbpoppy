"""Plot the change in DM distributions due to differing beampatterns."""
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, pprint

from tests.convenience import plot_aa_style, rel_path

BEAMPATTERNS = ['perfect', 'gaussian', 'airy-0', 'airy-4']
SIZE = 1e5


def get_data():
    """Get the population data."""
    # Construct population
    pop = CosmicPopulation(n_srcs=SIZE, n_days=1, name='standard_candle')
    pop.set_dist(model='sfr', z_max=2.5, H_0=67.74, W_m=0.3089, W_v=0.6911)
    pop.set_dm_host(model='constant', value=100)
    pop.set_dm_igm(model='ioka', slope=1000, std=None)
    pop.set_dm_mw(model='ne2001')
    pop.set_emission_range(low=10e6, high=10e9)
    pop.set_lum(model='constant', value=1e36)
    pop.set_w(model='constant', value=1.)
    pop.set_si(model='constant', value=0)
    pop.generate()

    # Survey population
    pops = {}
    for b in BEAMPATTERNS:
        pprint(f'Surveying with {b} beampattern')
        n_s = 0
        bp = b
        if b.startswith('airy'):
            bp, n_s = b.split('-')
            n_s = int(n_s)

        survey = Survey(name='perfect-small')
        # Prevent beam from getting larger than the sky
        survey.set_beam(model=bp, n_sidelobes=n_s, size=10)
        surv_pop = SurveyPopulation(pop, survey)
        print(surv_pop.source_rate)
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


if __name__ == '__main__':
    pops = get_data()
    plot_dm(pops)
