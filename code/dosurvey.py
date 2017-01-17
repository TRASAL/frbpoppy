from population import Population
from survey import Survey


def run(pop, survey_list):
    """Run any surveys and detect FRB sources"""

    # List of survey populations
    survey_pops = []

    for surv in survey_list:
        s = Survey(surv)

        survpop = Population()

        for src in pop.sources:

            snr = s.calc_snr(src, pop)
