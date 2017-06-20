import copy

from population import Population, unpickle
from survey import Survey


def observe(population,
            survey_name,
            pattern='gaussian',
            return_pop=True,
            scat=False,
            scint=False):
    """Run survey to detect FRB sources

    Args:
        population (Population): Population class of FRB sources to observe
        survey_name (str): Name of survey file with which to observe
        pattern (str): Gain pattern, either 'gaussian' (default) or 'airy'
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
    pop = unpickle(population.name)

    s = Survey(survey_name, pattern=pattern)
    surv_pop = Population()
    surv_pop.name = survey_name

    for src in pop.sources:

        # Detection flag
        det = False

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

            # Check if repeat FRBs are within an integration time
            if frb.time:
                if frb.time > s.t_obs:
                    s.frb_rates.out += 1
                    continue

            # Calculate observing properties such as the signal to noise ratio
            s.obs_prop(frb, src, pop)

            # Add scintillation
            if scint:
                s.scint(frb, src)

            # Check whether it has been detected
            if frb.snr > s.snr_limit:
                s.frb_rates.det += 1
                if not src.detected:
                    s.src_rates.det += 1
                    det = True
            else:
                s.frb_rates.faint += 1

            # If above 1 Jy
            if frb.s_peak > 1.0:
                s.frb_rates.jy += 1
                if not src.detected:
                    s.src_rates.jy += 1

            if det:
                src.detected = True

        if src.detected:
            surv_pop.add(src)
        else:
            s.src_rates.faint += 1

    # Scale rates according to length of survey etc
    s.scale_rates(pop)
    s.rates(pop)

    # Return population or survey
    if return_pop:
        return surv_pop
    else:
        return s
