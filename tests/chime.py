"""Determine whether frbpoppy can explain CHIME results."""
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np

from frbpoppy import Frbcat, RepeaterPopulation, Survey, SurveyPopulation
from frbpoppy import split_pop

from convenience import hist, plot_aa_style, rel_path


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
    r = RepeaterPopulation(1e5,
                           n_days=1,
                           dm_host_model='gaussian',
                           dm_host_mu=100,
                           dm_host_sigma=0,
                           dm_igm_index=1000,
                           dm_igm_sigma=0,
                           dm_mw_model='zero',
                           emission_range=[10e6, 10e9],
                           lum_range=[1e42, 1e45],
                           lum_index=0,
                           n_model='vol_co',
                           alpha=-1.5,
                           w_model='lognormal',
                           w_range=[1., 1.],
                           w_mu=0.1,
                           w_sigma=0.5,
                           si_mu=-1.4,
                           si_sigma=0.,
                           z_max=2.5,
                           lum_rep_model='independent',
                           lum_rep_sigma=1e3,
                           si_rep_model='same',
                           si_rep_sigma=0.1,
                           times_rep_model='even',
                           w_rep_model='independent',
                           w_rep_sigma=0.05,
                           generate=True)

    s = Survey('chime', strategy='follow-up', n_days=1)
    s.gain_pattern = 'perfect'
    s.snr_limit = 1.

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
