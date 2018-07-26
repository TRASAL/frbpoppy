"""Class to generate a survey population of FRBs."""
from copy import deepcopy

from frbpoppy.log import pprint
from frbpoppy.population import Population
from frbpoppy.survey import Survey
from frbpoppy.rates import NumberOf


class SurveyPopulation(Population):
    """Class to create a survey population of FRBs."""

    def __init__(self, cosmic_pop, survey, scat=False, scin=False):
        """
        Run a survey to detect FRB sources.

        Args:
            cosmic_pop (Population): Population class of FRB sources to observe
            survey (Survey): Survey class with which to observe
            scat (bool, optional): Whether to include scattering in signal to
                noise calculations.
            scin (bool, optional): Whether to apply scintillation to
                observations.
        """
        # Set up population
        Population.__init__(self)
        self.cosmic_pop = cosmic_pop
        self.survey = survey

        self.name = self.survey.name
        self.time = self.cosmic_pop.time
        self.vol_co_max = self.cosmic_pop.vol_co_max

        # TODO ADAPT
        # Copy population so that it can be observed multiple times
        # if not pop_path:
        #     pop = unpickle(population.name)
        # else:
        #     pop = unpickle(filename=pop_path)

        # Set rate counters
        self.num = NumberOf()

        for src in cosmic_pop.sources:

            # Check whether source is in region
            if not self.survey.in_region(src):
                self.num.srcs.out += 1
                self.num.frbs.out += src.n_frbs
                continue

            # Create a copy of the sourcce so it can be observed multiple times
            src = deepcopy(src)

            # Calculate dispersion measure across single channel, with error
            self.survey.dm_smear(src)

            # Set scattering timescale
            if scat:
                self.survey.scat(src)

            # Calculate total temperature
            self.survey.calc_Ts(src)

            for frb in src.frbs:

                # Check if repeat FRBs are within an integration time
                if frb.time:
                    if frb.time > self.survey.t_obs:
                        self.num.frbs.out += 1
                        continue

                # Calculate observing properties such as the S/R ratio
                self.survey.obs_prop(frb, src, pop)

                # Add scintillation
                if scin:
                    self.survey.scint(frb, src)

                # Check whether it has been detected
                if frb.snr > self.survey.snr_limit:
                    self.num.frbs.det += 1
                    if not src.detected:
                        self.num.srcs.det += 1
                        src.detected = True
                else:
                    self.num.frbs.faint += 1

            if src.detected:
                self.srcs.add(src)
            else:
                self.num.srcs.faint += 1

        # TODO ADAPT
        # s.rates(surv_pop, output=output)
