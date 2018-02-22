"""Run Monte Carlo over frbpoppy input parameters."""

from monte_carlo import MonteCarlo

mc = MonteCarlo()

mc.surveys = {'APERTIF': 'APERTIF',
              'PMSURV': 'parkes',
              'HTRU': 'parkes',
              'UTMOST-1D': 'UTMOST'}
            #   'PMSURV': 'parkes',
            #   'ASKAP-INCOH': 'ASKAP',
            #   'ASKAP-FLY': 'ASKAP',
            #   'GBT': 'GBT',
            #   'PALFA': 'arecibo',
            #   'ARECIBO-SPF': 'arecibo',
            #   'ALFABURST': 'arecibo'

# Set a limited run
mc.days = 60
mc.dm_host.limits(0, 100, 50, 100)
mc.dm_igm_slope.limits(1000, 1400, 100, 1200)
mc.freq_max.limits(1e7, 1e10, 1, 1e10, log=True)
mc.freq_min.limits(1e7, 1e10, 1, 1e7, log=True)
mc.lum_bol_max.limits(1e30, 1e60, 10, 1e50, log=True)
mc.lum_bol_min.limits(1e30, 1e60, 10, 1e40, log=True)
mc.lum_bol_slope.limits(-1.7, -1.3, 0.1, -1.5)
mc.n_day.limits(4000, 6000, 1000, 5000)
mc.rep.limits(0.0, 0.01, 0.01, 0)
mc.si_mean.limits(-1.6, -1.4, 0.1, -1.5)
mc.si_sigma.limits(0.0, 0.1, 0.1, 0.0)
mc.w_int_max.limits(0, 10, 0.5, 5)
mc.w_int_min.limits(0, 10, 0.5, 1)

mc.run()
