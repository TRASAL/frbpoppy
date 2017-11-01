"""Run Monte Carlo over frbpoppy input parameters."""

from monte_carlo import MonteCarlo
# from plot_mc import Plot

mc = MonteCarlo()
mc.dm_host.limits(0, 300, 100, 100)
mc.dm_host.limits(0, 200, 100, 100)
mc.dm_igm_slope.limits(1000, 1400, 200, 1200)
mc.freq_max.limits(10e5, 10e10, 1, 1e9, log=True)
mc.freq_min.limits(10e5, 10e10, 1, 1e6, log=True)
mc.lum_bol_max.limits(1e30, 1e60, 10, 1e50, log=True)
mc.lum_bol_min.limits(1e30, 1e60, 10, 1e40, log=True)
mc.lum_bol_slope.limits(0.5, 1.5, 0.5, 1.)
mc.n_day.limits(8000, 12000, 2000, 10000)
mc.rep.limits(0.0, 0.01, 0.01, 0)
mc.si_mean.limits(-2.0, -1.5, 0.1, -1.4)
mc.si_sigma.limits(0.0, 0.1, 0.1, 0.0)
mc.w_int_max.limits(0, 2, 1., 2)
mc.w_int_min.limits(0, 2, 1., 1)
mc.run()

# Plot().mc()
