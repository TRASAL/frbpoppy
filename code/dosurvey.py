from population import Population
from survey import Survey


def run(pop, survey_list):
    """Run any surveys and detect FRB sources"""

    # List of survey populations
    survey_pops = []

    for surv in survey_list:

        s = Survey(surv)
        surv_pop = Population()

        # Counters
        n_det = 0

        for src in pop.sources:

            # Calculate signal to noise ratio
            snr = s.calc_snr(src, pop)

            if snr > s.snr_limit:
                n_det += 1
                src.snr = snr_limit
                surv_pop.sources.append(src)
