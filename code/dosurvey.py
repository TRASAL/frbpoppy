from population import Population
from survey import Survey


def observe(pop, survey_name):
    """
    Run survey to detect FRB sources

    Args:
        pop (class): Population class of FRB sources to observe
        survey_name (str): Name of survey with which to observe
    Returns:
        surv_pop (class): Observed survey population
    """

    s = Survey(survey_name)
    surv_pop = Population()
    surv_pop.name = survey_name

    for src in pop.sources:

        # Calculate signal to noise ratio
        snr = s.calc_snr(src, pop)

        if snr > s.snr_limit:

            # Note that source has been detected
            s.n_det += 1
            src.snr = snr
            surv_pop.sources.append(src)

        elif snr == -2.0:
            s.n_out += 1
        else:
            s.n_faint += 1

    print(len(pop.sources), s.n_det, s.n_faint, s.n_out)

    return surv_pop
