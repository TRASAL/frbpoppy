"""Run Monte Carlo over frbpoppy input parameters."""

from monte_carlo import MonteCarlo
from plot_mc import Plot
from frbcat import get_frbcat

# mc = MonteCarlo()
# mc.dm_host.limits(0, 300, 50, 100)
# mc.dm_host.limits(0, 200, 50, 100)
# mc.dm_igm_slope.limits(1000, 1400, 100, 1200)
# mc.freq_max.limits(10e5, 10e10, 1, 10e9, log=True)
# mc.freq_min.limits(10e5, 10e10, 1, 10e6, log=True)
# mc.lum_bol_max.limits(1e30, 1e60, 10, 1e50, log=True)
# mc.lum_bol_min.limits(1e30, 1e60, 10, 1e40, log=True)
# mc.lum_bol_slope.limits(0.5, 1.5, 0.5, 1.)
# mc.n_day.limits(2000, 14000, 6000, 10000)
# mc.rep.limits(0.0, 0.1, 0.5, 0.05)
# mc.si_mean.limits(-2.0, -1, 0.5, -1.4)
# mc.si_sigma.limits(0.0, 0.5, 0.1, 0.0)
# mc.w_int_max.limits(0, 5, 1., 5)
# mc.w_int_min.limits(0, 5, 1., 1)
# mc.run()

Plot().mc()
