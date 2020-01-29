"""Determine whether frbpoppy can explain CHIME results."""
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np

from frbpoppy import Frbcat, CosmicPopulation, Survey, SurveyPopulation
from frbpoppy import split_pop

from convenience import hist, plot_aa_style, rel_path

N_DAYS = 100


def get_frbcat_data():
    """Get all chime data from frbcat.

    Returns:
        dict: Two keys 'r' for repeater and 'o' for one-offs. Each
            with entries for 'dm' and 'snr'

    """
    fc = Frbcat(frbpoppy=False, repeaters=True)
    chime_df = fc.df[fc.df.telescope == 'chime']
    chime_df = chime_df.sort_values(by=['frb_name'])

    frbcat = {'r': {}, 'o': {}}

    # Chime one-offs
    chime_o = chime_df.drop_duplicates(subset=['frb_name'], keep=False)

    # Chime repeaters
    chime_r = chime_df.loc[chime_df['frb_name'].duplicated(), :]

    # One DM value per repeater (used the average between bursts)
    frbcat['r']['dm'] = chime_r.groupby('frb_name').mean().reset_index().dm
    frbcat['o']['dm'] = chime_o.dm

    # All the different SNRs per repeater (or one_offs)
    frbcat['r']['snr'] = chime_r.snr
    frbcat['o']['snr'] = chime_o.snr

    # Number of repeaters vs one offs
    frbcat['r']['n'] = len(frbcat['r']['dm'])
    frbcat['o']['n'] = len(frbcat['o']['dm'])

    return frbcat


def get_frbpoppy_data():
    """Get frbpoppy data."""
    r = CosmicPopulation(1e5, n_days=N_DAYS, repeaters=True)
    r.set_dist(model='vol_co', z_max=2.5)
    r.set_dm_host(model='gauss', mu=100, sigma=0)
    r.set_dm_igm(model='ioka', slope=1000, sigma=0)
    r.set_dm(mw=False, igm=True, host=True)
    r.set_emission_range(low=100e6, high=10e9)
    r.set_lum(model='powerlaw', per_source='different',
              low=1e42, high=1e45, power=0)
    r.set_si(model='gauss', mu=-1.4, sigma=0)
    r.set_w(model='lognormal', per_source='different', mu=0.1, sigma=0.5)
    r.set_time(model='regular', lam=2)
    r.generate()

    s = Survey('chime', n_days=N_DAYS)
    s.set_beam(model='chime')

    surv_pop = SurveyPopulation(r, s)

    # Split population into seamingly one-off and repeater populations
    mask = ((~np.isnan(surv_pop.frbs.time)).sum(1) > 1)
    pop_ngt1, pop_nle1 = split_pop(surv_pop, mask)
    pop_ngt1.name += ' (> 1 burst)'
    pop_nle1.name += ' (1 burst)'

    frbpop = {'r': {}, 'o': {}}
    for i, pop in enumerate((pop_ngt1, pop_nle1)):
        t = 'o'
        if i == 0:
            t = 'r'
        frbpop[t]['dm'] = pop.frbs.dm
        frbpop[t]['snr'] = pop.frbs.snr

    return frbpop


def plot(frbcat, frbpop):
    """Plot distributions."""
    # Change working directory
    plot_aa_style(cols=2)

    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

    ax1.set_xlabel(r'DM ($\textrm{pc}\ \textrm{cm}^{-3}$)')
    ax1.set_ylabel('Fraction')
    ax2.set_xlabel(r'S/N')
    plt.xscale('log')

    # Set colours
    cmap = plt.get_cmap('tab10')([0, 1])

    # Plot dm distribution
    for i, p in enumerate((frbcat, frbpop)):
        for t in ['r', 'o']:

            # Line style
            linestyle = 'solid'
            label = 'one-offs'
            if t == 'r':
                linestyle = 'dashed'
                label = 'repeaters'

            bins = np.linspace(0, 3000, 10)
            ax1.step(*hist(p[t]['dm'], norm='max', bins=bins),
                     where='mid', linestyle=linestyle, label=label,
                     color=cmap[i])

            # Plot SNR distribution
            bins = np.logspace(-1, 6, 10)
            ax2.step(*hist(p[t]['snr'], norm='max', bins=bins),
                     where='mid', linestyle=linestyle, label=label,
                     color=cmap[i])

    # Set up layout options
    f.subplots_adjust(hspace=0)
    f.subplots_adjust(wspace=0.07)

    # Add legend elements
    elements = []

    def patch(color):
        return Patch(facecolor=color, edgecolor=color)

    elements.append((patch(cmap[0]), 'Frbcat'))
    elements.append((patch(cmap[1]), 'Frbpoppy'))
    # # Add gap in legend
    # elements.append((Line2D([0], [0], color='white'), ''))
    # Add line styles
    elements.append((Line2D([0], [0], color='gray'), 'One-offs'))
    elements.append((Line2D([0], [0], color='gray', linestyle='dashed'),
                     'Repeaters'))

    lines, labels = zip(*elements)

    lgd = plt.figlegend(lines, labels, loc='upper center', ncol=4,
                        framealpha=1,  bbox_to_anchor=(0.485, 1.04),
                        columnspacing=1.1, handletextpad=0.3)
    # plt.tight_layout()
    # plt.legend(lines, labels, bbox_to_anchor=(0.5, 1.04), loc="upper center")

    path = rel_path(f'./plots/frbcat_chime.pdf')
    plt.savefig(path, bbox_extra_artists=(lgd,), bbox_inches='tight')


def main():
    """Run code."""
    frbcat = get_frbcat_data()
    frbpop = get_frbpoppy_data()
    plot(frbcat, frbpop)


if __name__ == '__main__':
    main()
