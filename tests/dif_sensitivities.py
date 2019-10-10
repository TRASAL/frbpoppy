import numpy as np
import matplotlib.pyplot as plt
from frbpoppy import RepeaterPopulation, Survey, SurveyPopulation, plot
from frbpoppy import split_pop, pprint
from frbpoppy import galacticops as go

N = 5000


def limit_ra_dec(pop, pointings):
    """Doesn't work at boundary coordinates."""

    pprint('Limiting coordinates')

    def sample(n_gen):
        u = np.random.uniform
        ra =  u(0, 360, n_gen)
        dec = np.rad2deg(np.arccos(u(-1, 1, n_gen))) - 90
        return ra, dec

    def accept(ra, dec):
        coords = np.full(len(ra), False)
        r = np.sqrt(40/np.pi)
        for p_ra, p_dec in pointings:
            limits = go.separation(ra, dec, p_ra, p_dec) < r
            coords[limits] = True
        return coords

    # Limit population to smaller area
    # Sample RA, dec
    ra, dec = sample(pop.n_gen)
    mask = accept(ra, dec)
    reject, = np.where(~mask)
    while reject.size > 0:
        fill_ra, fill_dec = sample(reject.size)
        mask = accept(fill_ra, fill_dec)
        ra[reject[mask]] = fill_ra[mask]
        dec[reject[mask]] = fill_dec[mask]
        reject = reject[~mask]

    frbs = pop.frbs
    frbs.ra = ra
    frbs.dec = dec

    # Convert to galactic coordinates
    frbs.gl, frbs.gb = go.radec_to_lb(frbs.ra, frbs.dec, frac=True)
    pop.gen_gal_coords()

    return pop


def main():

    r = RepeaterPopulation(N,
                           days=1,
                           dm_host_model='gaussian',
                           dm_host_mu=100,
                           dm_host_sigma=0,
                           dm_igm_index=1000,
                           dm_igm_sigma=0,
                           dm_mw_model='zero',
                           emission_range=[10e6, 10e9],
                           lum_range=[1e38, 1e38],
                           lum_index=0,
                           n_model='vol_co',
                           alpha=-1.5,
                           w_model='uniform',
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

    s = Survey('perfect')
    s.gain_pattern = 'perfect'

    # Setup pointings
    n_p = 1  # # of pointings
    decs = np.linspace(s.dec_min, s.dec_max, n_p+2)[1:n_p+1]
    ras = np.linspace(s.ra_min, s.ra_max, n_p+2)[1:n_p+1]
    s.pointings = list(zip(ras, decs))
    pop = limit_ra_dec(r, s.pointings)
    pop.name = 'Cosmic Population'


    s.snr_limit = 1
    surv_pop_low_snr = SurveyPopulation(r, s)
    surv_pop_low_snr.name = 'Low SNR'

    s.snr_limit = 10
    surv_pop_high_snr = SurveyPopulation(r, s)
    surv_pop_high_snr.name = 'High SNR'

    pops = [pop]
    # Split population into seamingly one-off and repeater populations
    for pop in (surv_pop_low_snr, surv_pop_high_snr):
        mask = ((~np.isnan(pop.frbs.time)).sum(1)>1)
        pop_ngt1, pop_nle1 = split_pop(pop, mask)
        pop_ngt1.name += ' (> 1 burst)'
        pop_nle1.name += ' (1 burst)'
        pops.append(pop_nle1)
        pops.append(pop_ngt1)


    plot(*pops, frbcat=False, mute=False, show=True)


if __name__ == '__main__':
    main()
