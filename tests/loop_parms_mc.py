"""Test the looping over all parameters set in a MC run."""
from frbpoppy.monte_carlo import MonteCarlo

mc = MonteCarlo()

# Set up a small run
mc.surveys = {'HTRU': 'parkes',
              'UTMOST-1D': 'UTMOST'}

# Set a limited run
mc.days = 30
mc.dm_host.limits(0, 100, 100, 100)
mc.dm_igm_index.limits(1000, 1200, 100, 1200)
mc.freq_max.limits(1e7, 1e10, 3, 1e10, log=True)
mc.freq_min.limits(1e7, 1e10, 3, 1e7, log=True)
mc.lum_bol_max.limits(1e30, 1e60, 10, 1e40, log=True)
mc.lum_bol_min.limits(1e30, 1e60, 10, 1e30, log=True)
mc.lum_bol_slope.limits(-2.0, -1.5, 0.5, -1.5)
mc.n_day.limits(4000, 6000, 1000, 5000)
mc.rep.limits(0.0, 0.01, 0.01, 0)
mc.si_mu.limits(-1.6, -1.4, 0.1, -1.5)
mc.si_sigma.limits(0.0, 0.1, 0.1, 0.0)
mc.w_int_max.limits(0, 2, 1, 2.0)
mc.w_int_min.limits(0, 2, 1, 1.0)

mc.extension = 'test_small'

pp = mc.possible_pars()

print(pp)
