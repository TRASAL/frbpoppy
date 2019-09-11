import numpy as np
from frbpoppy import RepeaterPopulation, Survey, SurveyPopulation
from frbpoppy import split_pop, pprint
from frbpoppy import galacticops as go

from askap_repeaters import limit_ra_dec

N = 10000


def shift_nan_to_right(arr):
    n = arr.shape[1]
    mask = np.isnan(arr)
    idx = mask.argmax(1)
    idx[~mask.any(1)] = n
    arr[idx[:,None] <= np.arange(n)] = np.nan
    return arr

def plot_time_dif(pops):
    """Plot the difference in time between various FRBs."""
    import matplotlib.pyplot as plt

    f, (ax1, ax2) = plt.subplots(1, 2)

    for pop in pops:
        time = pop.frbs.time.flatten()
        time = time[~np.isnan(time)]
        weights = np.ones_like(time)/len(time)
        ax1.hist(time, label=pop.name, weights=weights,
                 histtype='step', bins=20, alpha=0.5, linewidth=2)

        time = shift_nan_to_right(pop.frbs.time)
        t_diff = np.diff(time).flatten()
        t_diff = t_diff[~np.isnan(t_diff)]
        weights = np.ones_like(t_diff)/len(t_diff)
        if t_diff.any():
            ax2.hist(t_diff, weights=weights, label=pop.name,
                     histtype='step', bins=20, alpha=0.5,
                     linewidth=2)


    ax1.set_xlabel(r'$t_{bursts}$')
    ax1.set_ylabel('Fraction of bursts')
    ax2.set_xlabel(r'$\Delta t_{bursts}$')
    ax2.set_yscale('log')
    ax2.legend()
    plt.tight_layout()
    plt.savefig('./plots/times.pdf')


def main():

    r = RepeaterPopulation(N,
                           days=1,
                           dm_host_model='normal',
                           dm_host_mu=100,
                           dm_host_sigma=0,
                           dm_igm_index=1000,
                           dm_igm_sigma=0,
                           dm_mw_model='zero',
                           emission_range=[10e6, 10e9],
                           lum_range=[1e43, 1e43],
                           lum_index=0,
                           n_model='vol_co',
                           alpha=-1.5,
                           w_model='uniform',
                           w_range=[1., 1.],
                           w_mu=0.1,
                           w_sigma=0.5,
                           si_mu=-1.4,
                           si_sigma=0.,
                           z_max=0.1,
                           lum_rep_model='same',
                           lum_rep_sigma=1e3,
                           si_rep_model='same',
                           si_rep_sigma=0.1,
                           times_rep_model='even',
                           w_rep_model='same',
                           w_rep_sigma=0.05,
                           generate=True)

    s = Survey('askap-incoh')
    s.gain_pattern = 'perfect'
    s.snr_limit = 0.1

    # Setup pointings
    decs = np.linspace(s.dec_min, s.dec_max, 5)[1:4]
    ras = np.linspace(s.ra_min, s.ra_max, 5)[1:4]
    s.pointings = list(zip(ras, decs))
    pop = limit_ra_dec(r, s.pointings)
    pop.name = 'Cosmic'

    surv_pop = SurveyPopulation(r, s)
    surv_pop.name = 'Askap'

    plot_time_dif([pop, surv_pop])

if __name__ == '__main__':
    main()
