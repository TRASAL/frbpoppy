"""Test whether distributions for simple populations hold true."""
import os
import sys
sys.path.append("..")
from do_plot import plot
from do_populate import generate
from do_survey import observe

folder = os.path.dirname(os.path.realpath(__file__)) + '/'

days = 1
n_per_day = 10000

# Generate FRB population
population = generate(n_per_day*days,
                      days=days,
                      lum_dist_pars=[1e42, 1e42, -1.0],
                      z_max=0.1,
                      dm_pars=[0, 1200],
                      electron_model='zero',
                      emission_pars=[10e6, 10e9],
                      pulse=[5, 5],
                      si_pars=[0., 0.],
                      repeat=0.0,
                      save_path=folder + 'simple_pop.p')

population.save(out=folder + 'pop_initial.csv')

# Observe FRB populations
surveys = ['APERTIF']#,
        #    'HTRU',
        #    'UTMOST-1D']

pop_paths = []

for s in surveys:
    surv_pop = observe('_', s, pop_path=folder + 'simple_pop.p')
    pop_path = folder + 'pop_' + s.lower() + '.csv'
    surv_pop.save(out=pop_path)
    pop_paths.append(pop_path)

pop_paths = [folder + 'pop_' + 'initial.csv', *pop_paths]

# Plot populations
plot(files=pop_paths, mute=False)
