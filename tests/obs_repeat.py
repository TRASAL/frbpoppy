"""Test FRB population consisting of repeaters."""
from frbpoppy.do_populate import generate
from frbpoppy.do_plot import plot
from frbpoppy.do_survey import observe
from frbpoppy.paths import paths

MAKE = False

if MAKE:
    days = 7
    n_per_day = 5000

    # Generate population with full repeaters
    pop_repeat = generate(n_per_day*days,
                           days=days,
                           repeat=1.,
                           name='repeat')

    surv_pop_repeat = observe(pop_repeat, 'APERTIF', gain_pattern='apertif')
    surv_pop_repeat.name = 'obs-repeat'

    plot_args = [pop_repeat, surv_pop_repeat]

    # Plot populations
    plot(*plot_args)

else:
    pops = ['repeat', 'obs-repeat']
    plot_args = [f'{paths.populations()}population_{n}.csv' for n in pops]

    # Plot populations
    plot(files=plot_args)
