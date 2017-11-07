"""Run Monte Carlo over frbpoppy input parameters."""

# from monte_carlo import MonteCarlo
#
# mc = MonteCarlo()
#
# # Set a limited run
# mc.days = 30
# mc.dm_host.limits(0, 100, 10, 100)
# mc.dm_igm_slope.limits(1000, 1400, 50, 1200)
# mc.freq_max.limits(10e5, 10e10, 1, 10e10, log=True)
# mc.freq_min.limits(10e5, 10e10, 1, 10e5, log=True)
# mc.lum_bol_max.limits(1e30, 1e60, 10, 1e60, log=True)
# mc.lum_bol_min.limits(1e10, 1e30, 10, 1e20, log=True)
# mc.lum_bol_slope.limits(-0.9, -1.3, 0.1, -1.1)
# mc.n_day.limits(4000, 6000, 1000, 5000)
# mc.rep.limits(0.0, 0.05, 0.01, 0)
# mc.si_mean.limits(-2.0, -1.5, 0.1, -1.6)
# mc.si_sigma.limits(0.0, 0.5, 0.1, 0.0)
# mc.w_int_max.limits(0, 2, 0.2, 2)
# mc.w_int_min.limits(0, 2, 0.2, 1)
#
# mc.run()

import frbcat

df = frbcat.get_frbcat()
print(df.info())
