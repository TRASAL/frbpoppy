"""Create a standard CosmicPopulation to observe."""

def gen_standard():
    """Generate a standard population."""
    from frbpoppy import CosmicPopulation

    n_gen = int(4e5)

    pop = CosmicPopulation(n_gen,
                         days=1,
                         name='standard',
                         H_0=69.6,
                         W_m=0.286,
                         W_v=0.714,
                         dm_host_model='normal',
                         dm_host_mu=100,
                         dm_host_sigma=0,
                         dm_igm_index=1000,
                         dm_igm_sigma=None,
                         dm_mw_model='ne2001',
                         emission_range=[10e6, 10e9],
                         lum_range=[1e40, 1e45],
                         lum_index=0.,
                         n_model='vol_co',
                         alpha=-1.5,
                         pulse_model='uniform',
                         pulse_range=[1., 1.],
                         pulse_mu=1.6,
                         pulse_sigma=1.,
                         si_mu=-1.4,
                         si_sigma=1.,
                         z_max=2.5)

    pop.save()

    return pop
