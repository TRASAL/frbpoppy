import copy

from population import Population
from survey import Survey


def observe(population,
            survey_name,
            return_pop=True,
            scat=False,
            scint=False):
    """Run survey to detect FRB sources

    Args:
        population (Population): Population class of FRB sources to observe
        survey_name (str): Name of survey file with which to observe
        return_pop (bool, optional): Whether to return a population or survey
            class. Primarily intended for debugging. Defaults to True
        scat (bool, optional): Whether to include scattering in signal to noise
            calculations. Defaults to False
        scint (bool, optional): Whether to apply scintillation to observations.
            Defaults to False

    Returns:
        surv_pop (Population): Observed survey population
    """

    # Copy population so that it can be observed multiple times
    pop = copy.deepcopy(population)

    s = Survey(survey_name)
    surv_pop = Population()
    surv_pop.name = survey_name

    for src in pop.sources:

        # Check whether source is in region
        if not s.in_region(src):
            s.src_rates.out += 1
            s.frb_rates.out += src.n_frbs
            continue

        # Calculate dispersion measure across single channel, with error
        s.dm_smear(src)

        # Set scattering timescale
        if scat:
            s.scat(src)

        # Calculate total temperature
        s.calc_Ts(src)

        for frb in src.frbs:

            # Calculate observing properties such as the signal to noise ratio
            s.obs_prop(frb, src, pop)

            # Add scintillation
            if scint:
                s.scint(frb, src)

            if frb.snr > s.snr_limit:
                # Note that frb has been detected
                s.frb_rates.det += 1

                if not src.detected:
                    s.src_rates.det += 1
                    src.detected = True

            else:
                s.frb_rates.faint += 1

        if src.detected:
            surv_pop.add(src)
        else:
            s.src_rates.faint += 1

    s.rates()

    # Return population or survey
    if return_pop:
        return surv_pop
    else:
        return s
