"""Plot the DM distribution obtained with frbpoppy against frbcat results."""
import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, unpickle
from frbpoppy import Frbcat, pprint

from gen_standard import gen_standard
from plot_frbpoppy_parkes import hist, plot_dists

MAKE = False
OBSERVE = False

def main():

    if MAKE:
        pop = gen_standard()

    if OBSERVE:

        if not MAKE:
            pop = unpickle('standard')

        survey = Survey('askap-fly', gain_pattern='airy')
        surv_pop = SurveyPopulation(pop, survey)
        surv_pop.name = 'standard_askap'
        surv_pop.save()

    else:
        surv_pop = unpickle('standard_askap')

    plot_dists(surv_pop, 'askap')

if __name__ == '__main__':
    main()
