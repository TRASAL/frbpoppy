from population import Population
from survey import Survey


def observe(pop,
            survey_name,
            return_pop=True,
            scat=False,
            scint=False,):
    """Run survey to detect FRB sources

    Args:
        pop (Population): Population class of FRB sources to observe
        survey_name (str): Name of survey file with which to observe
        return_pop (bool): Whether to return a population or survey class.
            Primarily intended for debugging. Defaults to True
        scat (bool): Whether to include scattering in signal to noise
            calculations. Defaults to False
        scint (bool): Whether to apply scintillation to observations. Defaults
            to False

    Returns:
        surv_pop (Population): Observed survey population
    """

    s = Survey(survey_name)
    surv_pop = Population()
    surv_pop.name = survey_name

    for src in pop.sources:

        # Calculate observing properties such as the signal to noise ratio,
        # effective pulse width etc.
        snr, w_eff, s_peak, fluence = s.obs_prop(src, pop, scat=scat)

        # Check whether source is outside survey region
        if snr == -2.0:
            s.n_out += 1
            continue

        # Add scintillation
        if scint:
            snr = s.scint(src, snr)

        if snr > s.snr_limit:
            # Note that source has been detected
            s.n_det += 1
            src.detected = True
            src.snr = snr
            src.w_eff = w_eff
            src.s_peak = s_peak
            src.fluence = fluence
            surv_pop.sources.append(src)
        else:
            s.n_faint += 1

    s.result()

    if return_pop:
        return surv_pop
    else:
        return s
