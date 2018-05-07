"""Allow a survey to be run over a population of FRBs."""
from frbpoppy.population import Population, unpickle
from frbpoppy.survey import Survey
from frbpoppy.paths import paths


def observe(population,
            survey_name,
            gain_pattern='gaussian',
            sidelobes=1,
            equal_area=False,
            output=True,
            return_pop=True,
            return_survey=False,
            scat=False,
            scint=False,
            pop_path=''):
    """
    Run survey to detect FRB sources.

    Args:
        population (Population): Population class of FRB sources to observe
        survey_name (str): Name of survey file with which to observe
        gain_pattern (str, optional): Gain pattern, either 'gaussian' (default)
            'airy', 'tophat' or 'perfect'
        output (bool, optional): Whether to print the rates or not
        return_pop (bool, optional): Whether to return a population. Primarily
            intended for debugging. Defaults to True
        return_survey (bool, optional): Whether to return a survey. Primarily
            intended for debugging. Defaults to False
        scat (bool, optional): Whether to include scattering in signal to noise
            calculations. Defaults to False
        scint (bool, optional): Whether to apply scintillation to observations.
            Defaults to False
        pop_path (str): Give filename for a pickled population rather than
            using a population class

    Returns:
        surv_pop (Population): Observed survey population

    """
    # Copy population so that it can be observed multiple times
    if not pop_path:
        pop = unpickle(population.name)
    else:
        pop = unpickle(filename=pop_path)

    s = Survey(survey_name,
               gain_pattern=gain_pattern,
               sidelobes=sidelobes,
               equal_area=equal_area)
    surv_pop = Population()
    surv_pop.name = survey_name
    surv_pop.time = pop.time
    surv_pop.vol_co_max = pop.vol_co_max

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
    s.scale_rates(surv_pop)
    s.rates(surv_pop, output=output)

    # Return population or survey classes
    if return_pop and not return_survey:
        return surv_pop
    if return_pop and return_survey:
        return s, surv_pop
    if return_survey and not return_pop:
        return s
