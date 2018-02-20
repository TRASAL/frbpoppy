from do_populate import generate
from do_survey import observe
from do_plot import plot

# Generate FRB population
population = generate(30000,
                      days=3,
                      lum_dist_pars=[1e41, 1e45, -1.4],
                      z_max=2.5,
                      pulse=[0.1, 10],
                      repeat=0.05)

# Observe FRB populations
surveys = ['WHOLESKY',
           'APERTIF',
           'APERTIF-IAB-1',
           'APERTIF-IAB-10',
           'APERTIF-IAB-12',
           'PMSURV',
           'HTRU',
           'ASKAP-INCOH',
           'ASKAP-FLY',
           'GBT',
           'PALFA',
           'ARECIBO-SPF',
           'ALFABURST',
           'UTMOST-1D']

# Other options
# 'ASKAP-FLY',
# 'VLA-L-BAND',
# 'UTMOST-2D',
# 'MWA',
# 'LWA',
# 'MEERKAT',
# 'CHIME',
# 'DSA10',
# 'EMERLIN',
# 'VFASTR',

results = []

for s in surveys:
    results.append(observe(population, s))
# wholesky = observe(population, 'WHOLESKY')
# test = observe(population, 'TEST')

# Plot populations
# plot(population, *results, mute=False)
