"""Class to generate a survey population of FRBs."""
from copy import deepcopy
import math
import numpy as np

from frbpoppy.log import pprint
from frbpoppy.population import Population
from frbpoppy.repeater_pop import RepeaterPopulation
from frbpoppy.rates import Rates, scale


class SurveyPopulation(Population):
    """Class to create a survey population of FRBs."""

    def __init__(self, cosmic_pop, survey, scat=False, scin=False,
                 rate_limit=True):
        """
        Run a survey to detect FRB sources.

        Args:
            cosmic_pop (Population): Population class of FRB sources to observe
            survey (Survey): Survey class with which to observe
            scat (bool, optional): Whether to include scattering in signal to
                noise calculations.
            scin (bool, optional): Whether to apply scintillation to
                observations.
            rate_limit (bool, optional): Whether to limit detections by 1/(1+z)
                due to limitation in observing time
        """
        pprint(f'Surveying {cosmic_pop.name} with {survey.name}')
        # Stops RuntimeWarnings about nan values
        np.warnings.filterwarnings('ignore')

        # Set up population
        Population.__init__(self)

        # Set attributes
        self.name = survey.name
        self.time = cosmic_pop.time
        self.vol_co_max = cosmic_pop.vol_co_max
        self.frbs = deepcopy(cosmic_pop.frbs)
        self.rate = Rates()

        frbs = self.frbs

        # Check whether source is in region
        region_mask = survey.in_region(frbs)
        frbs.apply(region_mask)
        self.rate.out = np.size(region_mask) - np.count_nonzero(region_mask)

        # Calculate dispersion measure across single channel
        frbs.t_dm = survey.dm_smear(frbs)

        # Set scattering timescale
        if scat:
            frbs.t_scat = survey.calc_scat(frbs.dm)

        # Calculate total temperature
        frbs.T_sky, frbs.T_sys = survey.calc_Ts(frbs)

        # Calculate effective pulse width
        frbs.w_eff = survey.calc_w_eff(frbs)

        # Calculate peak flux density
        f_min = cosmic_pop.f_min
        f_max = cosmic_pop.f_max
        frbs.s_peak = survey.calc_s_peak(frbs, f_low=f_min, f_high=f_max)

        # Account for beam offset
        int_pro, offset = survey.intensity_profile(shape=frbs.s_peak.shape)
        if int_pro.ndim == frbs.s_peak.ndim:
            frbs.s_peak *= int_pro
        else:
            frbs.s_peak *= int_pro[:, None]
        frbs.offset = offset / 60.  # In degrees

        # Calculate fluence [Jy*ms]
        frbs.fluence = survey.calc_fluence(frbs)

        # Caculate Signal to Noise Ratio
        frbs.snr = survey.calc_snr(frbs)

        # Add scintillation
        if scin:
            # Not sure whether this can cope with 2D arrays
            frbs.snr = survey.calc_scint(frbs)

        # Check whether frbs would be above detection threshold
        snr_mask = (frbs.snr >= survey.snr_limit)
        frbs.apply(snr_mask)
        self.rate.faint = np.size(snr_mask) - np.count_nonzero(snr_mask)

        if isinstance(cosmic_pop, RepeaterPopulation):
            # Add a time filter on repeaters
            obs_strategy_mask = survey.in_observation(frbs)
            frbs.apply(obs_strategy_mask)
            frbs.clean_up()
        elif rate_limit is True:
            limit = 1/(1+frbs.z)
            rate_mask = np.random.random(len(frbs.z)) <= limit
            frbs.apply(rate_mask)
            self.rate.late = np.size(rate_mask) - np.count_nonzero(rate_mask)

        self.rate.det = len(frbs.snr)

        if not isinstance(cosmic_pop, RepeaterPopulation):
            # Scale rates by beam size etc if one-off population
            self.calc_rates(survey)

    def calc_rates(self, survey):
        """Calculate the relative detection rates."""
        # Calculate scaling factors for rates
        area_sky = 4*math.pi*(180/math.pi)**2   # In sq. degrees
        f_area = (survey.beam_size * self.rate.tot())
        inside = self.rate.det+self.rate.late+self.rate.faint
        f_area /= (inside*area_sky)
        f_time = 86400 / self.time

        # Saving scaling factors
        self.rate.days = self.time/86400  # seconds -> days
        self.rate.f_area = f_area
        self.rate.f_time = f_time
        self.rate.name = self.name
        self.rate.vol = self.rate.tot()
        self.rate.vol /= self.vol_co_max * (365.25*86400/self.time)

    def rates(self, scale_area=True, scale_time=False):
        """Scale the rates by the beam area and time."""
        r = self.rate
        if scale_area:
            r = scale(self.rate, area=True)
        if scale_time:
            r = scale(self.rate, time=True)
        return r

    def calc_logn_logs(self, parameter='fluence', min_p=None, max_p=None):
        """TODO. Currently unfinished."""
        parms = getattr(self.frbs, parameter)

        if min_p is None:
            f_0 = min(parms)
        else:
            f_0 = min_p
            parms = parms[parms >= min_p]

        if max_p is not None:
            parms = parms[parms <= max_p]

        n = len(parms)
        alpha = -1/((1/n)*sum([math.log(f/f_0) for f in parms]))
        alpha *= (n-1)/n  # Removing bias in alpha
        alpha_err = n*alpha/((n-1)*(n-2)**0.5)
        norm = n / (f_0**alpha)  # Normalisation at lowest parameter

        return alpha, alpha_err, norm


if __name__ == '__main__':
    from frbpoppy import CosmicPopulation, Survey
    cosmic = CosmicPopulation(10000, days=4)
    survey = Survey('apertif', gain_pattern='apertif')
    surv_pop = SurveyPopulation(cosmic, survey)
    print(surv_pop.rates())
