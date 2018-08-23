"""Class to generate a survey population of FRBs."""
from copy import deepcopy
import math
import numpy as np
from lmfit import Model

from frbpoppy.log import pprint
from frbpoppy.population import Population
from frbpoppy.rates import Rates, scale


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

        # Set rate counters
        self.frb_rates = Rates(rate_type='FRBs')
        self.src_rates = Rates(rate_type='Sources')

        for src in cosmic_pop.sources:

            # Check whether source is in region
            if not self.survey.in_region(src):
                self.src_rates.out += 1
                self.frb_rates.out += src.n_frbs
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
                        self.frb_rates.out += 1
                        continue

                # Calculate observing properties such as the S/R ratio
                args = (frb, src, self.cosmic_pop.f_min, self.cosmic_pop.f_max)
                self.survey.obs_prop(*args)

                # Add scintillation
                if scin:
                    self.survey.scint(frb, src)

                # Check whether it has been detected
                if frb.snr > self.survey.snr_limit:
                    self.frb_rates.det += 1
                    if not src.detected:
                        self.src_rates.det += 1
                        src.detected = True
                else:
                    self.frb_rates.faint += 1

            if src.detected:
                self.add(src)
            else:
                self.src_rates.faint += 1

        # Calculate scaling factors
        area_sky = 4*math.pi*(180/math.pi)**2   # In sq. degrees
        f_area = (self.survey.beam_size * self.src_rates.tot())
        f_area /= ((self.src_rates.det + self.src_rates.faint)*area_sky)
        f_time = 86400 / self.time

        # Saving scaling factors
        for r in (self.frb_rates, self.src_rates):
            r.days = self.time/86400
            r.f_area = f_area
            r.f_time = f_time
            r.name = self.name
            r.vol = r.tot() / self.vol_co_max * (365.25*86400/self.time)

    def rates(self, scale_area=True, scale_time=False, type='frbs'):
        """Adapt frb or source rates as needed."""
        if type in ('frb', 'frbs'):
            r = self.frb_rates
        elif type in ('src', 'srcs', 'sources'):
            r = self.src_rates

        if scale_area:
            r = scale(r, area=True)
        if scale_time:
            r = scale(r, time=True)

        return r

    def calc_logn_logs(self):
        """TODO. Currently unfinished."""
        f_lim = self.survey.calc_fluence_limit()
        f_lim = 5e-1
        fluences = [f for f in self.get('fluence') if f >= f_lim]

        # Bin up
        min_f = np.log10(min(fluences))
        max_f = np.log10(max(fluences))
        log_bins = np.logspace(min_f, max_f, 50)
        hist, edges = np.histogram(fluences, bins=log_bins)
        cum_hist = [sum(hist[i:]) for i in range(len(hist))]

        log_xs = ((np.log10(edges[:-1]) + np.log10(edges[1:])) / 2)
        xs = 10**log_xs

        def powerlaw(x, norm=1, alpha=-1):
            return norm*x**(alpha)

        def linear(x, norm=1, alpha=-1):
            """log(y) = log(norm) + alpha * log(x)"""
            return norm + alpha*x

        in_logspace = False

        if in_logspace:
            model = Model(linear)
            params = model.make_params()

            params['alpha'].min = -5.0
            params['alpha'].max = 0.5

            elems = zip(log_xs, np.log10(cum_hist))
            elems = [e for e in elems if e[1] > 0]
            log_xs, log_cum = zip(*elems)
            result = model.fit(log_cum, params, x=log_xs)

            alpha = result.params['alpha'].value
            alpha_err = result.params['alpha'].stderr
            norm = 10**(result.params['norm'].value)
        else:
            model = Model(powerlaw)
            params = model.make_params()

            params['alpha'].min = -5.0
            params['alpha'].max = 0.5

            result = model.fit(cum_hist, params, x=xs)

            alpha = result.params['alpha'].value
            alpha_err = result.params['alpha'].stderr
            norm = result.params['norm'].value

        return alpha, alpha_err, norm