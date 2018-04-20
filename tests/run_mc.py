"""Run Monte Carlo over frbpoppy input parameters."""

from frbpoppy.monte_carlo import MonteCarlo

mc = MonteCarlo()

mc.surveys = {'HTRU': 'parkes',
              'UTMOST-1D': 'UTMOST'}

# Set a limited run
mc.days = 90
mc.dm_host.limits(0, 100, 10, 100)
mc.dm_igm_slope.limits(600, 1400, 50, 1200)
mc.freq_max.limits(1e7, 1e10, 1, 1e10, log=True)
mc.freq_min.limits(1e7, 1e10, 1, 1e7, log=True)
mc.lum_bol_max.limits(1e30, 1e60, 5, 1e40, log=True)
mc.lum_bol_min.limits(1e30, 1e60, 5, 1e30, log=True)
mc.lum_bol_slope.limits(-2.0, -1.3, 0.1, -1.5)
mc.n_day.limits(4000, 6000, 1000, 5000)
mc.rep.limits(0.0, 0.01, 0.01, 0)
mc.si_mean.limits(-1.6, -1.4, 0.1, -1.5)
mc.si_sigma.limits(0.0, 0.1, 0.1, 0.0)
mc.w_int_max.limits(0, 2, 0.1, 1.0)
mc.w_int_min.limits(0, 2, 0.1, 0.1)

mc.run()
