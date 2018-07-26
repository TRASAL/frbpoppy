"""Test FRB sources following star forming rate."""
from frbpoppy.do_populate import generate
from frbpoppy.do_plot import plot
from frbpoppy.do_survey import observe

MAKE = True

days = 7
n_per_day = 5000

# Generate population following a uniform pulse distribution
pop_uniform = generate(n_per_day*days,
                       days=days,
                       pulse_model='uniform',
                       name='uniform')

# Generate population following log normal pulse distribution
pop_lognormal = generate(n_per_day*days,
                         days=days,
                         pulse_model='lognormal',
                         name='lognormal')

surv_pop_uniform = observe(pop_uniform, 'APERTIF', gain_pattern='apertif')
surv_pop_uniform.name = 'obs-uniform'

surv_pop_log = observe(pop_lognormal, 'APERTIF', gain_pattern='apertif')
surv_pop_log.name = 'obs-lognormal'

# Plot populations
plot(pop_uniform, pop_lognormal, surv_pop_uniform, surv_pop_log)
