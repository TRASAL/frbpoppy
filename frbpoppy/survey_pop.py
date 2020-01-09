"""Class to generate a survey population of FRBs."""
from copy import deepcopy
from tqdm import tqdm
import math
import numpy as np

from frbpoppy.log import pprint
from frbpoppy.population import Population
from frbpoppy.rates import Rates, scale
from frbpoppy.repeater_pop import RepeaterPopulation
import frbpoppy.galacticops as go


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
        self.scat = scat
        self.scin = scin
        self.survey = survey

        # Calculations differ for repeaters
        self.repeaters = False
        if isinstance(cosmic_pop, RepeaterPopulation):
            self.repeaters = True

            if scin is True:
                m = 'Scintillation is currently not implemented for '
                m += 'RepeaterPopulations'
                raise ValueError(m)

        frbs = self.frbs

        # Check whether source is in region
        region_mask = survey.in_region(frbs.ra, frbs.dec, frbs.gl, frbs.gb)
        frbs.apply(region_mask)
        self.rate.out = np.size(region_mask) - np.count_nonzero(region_mask)

        # Calculate dispersion measure across single channel
        frbs.t_dm = survey.dm_smear(frbs.dm)

        # Set scattering timescale
        if scat:
            frbs.t_scat = survey.calc_scat(frbs.dm)

        # Calculate total temperature
        frbs.T_sky, frbs.T_sys = survey.calc_Ts(frbs.gl, frbs.gb)

        # Calculate effective pulse width
        frbs.w_eff = survey.calc_w_eff(frbs.w_arr, frbs.t_dm, frbs.t_scat)

        # Calculate peak flux density
        frbs.s_peak = survey.calc_s_peak(frbs.si,
                                         frbs.lum_bol,
                                         frbs.z,
                                         frbs.dist_co,
                                         frbs.w_arr,
                                         frbs.w_eff,
                                         f_low=cosmic_pop.f_min,
                                         f_high=cosmic_pop.f_max)

        # Calculations differ whether dealing with repeaters or not
        if self.repeaters:
            self.det_repeaters()
        else:
            self.det_oneoffs()

        # Prevent additional memory usage
        self.survey = None

    def det_oneoffs(self):
        """Detect one-off frbs."""
        frbs = self.frbs
        survey = self.survey

        # Account for beam offset
        int_pro, offset = survey.calc_int_pro(shape=frbs.s_peak.shape)
        frbs.s_peak *= int_pro
        frbs.offset = offset  # [deg]

        # Calculate fluence [Jy*ms]
        frbs.fluence = survey.calc_fluence(frbs.s_peak, frbs.w_eff)

        # Calculate Signal to Noise Ratio
        frbs.snr = survey.calc_snr(frbs.s_peak, frbs.w_arr, frbs.T_sys)

        # Add scintillation
        if self.scin:

            # Ensure scattering has been calculated
            if not isinstance(frbs.t_scat, np.ndarray):
                frbs.t_scat = survey.calc_scat(frbs.dm)

            # Calculate signal to noise ratio after scattering
            frbs.snr = survey.calc_scint(frbs.t_scat, frbs.dist_co, frbs.gl,
                                         frbs.gb, frbs.snr)

        # Check whether frbs would be above detection threshold
        snr_mask = (frbs.snr >= survey.snr_limit)
        frbs.apply(snr_mask)
        self.rate.faint = np.size(snr_mask) - np.count_nonzero(snr_mask)

        # Distant frbs are redshifted out of your observing time
        limit = 1/(1+frbs.z)
        rate_mask = np.random.random(len(frbs.z)) <= limit
        frbs.apply(rate_mask)
        self.rate.late = np.size(rate_mask) - np.count_nonzero(rate_mask)

        self.rate.det = len(frbs.snr)

        # Calculate detection rates
        self.calc_rates(survey)

    def det_repeaters(self):
        """Detect repeating frbs."""
        frbs = self.frbs
        survey = self.survey

        # Set up a tuple of pointings if not given
        if survey.pointings is None:
            survey.gen_pointings()

        # t_obs in fractional days
        t_obs = survey.t_obs / 86400

        # Array with times of each pointing
        max_t = survey.n_days
        times = np.linspace(0, max_t, max_t/t_obs + 1)  # [days]
        lsts = times*360*(24/23.9344696) % 360  # Local sidereal time [deg]
        lsts += np.random.uniform(0, 360)  # Add random offset

        # Only keep bursts within survey time
        time_mask = (frbs.time <= times[-1])
        frbs.apply(time_mask)

        # Prepare for iterating over time
        p_i = 0  # Pointing index
        r = survey.calc_beam_radius()  # Beam pattern radius
        max_n_pointings = len(times[:-1])

        # Initialize some necessary arrays
        if frbs.w_eff.ndim == 2 or frbs.lum_bol.ndim == 2:
            sim_shape = frbs.time  # 2D
        else:
            sim_shape = frbs.lum_bol  # 1D

        frbs.fluence = np.full_like(sim_shape, np.nan)
        frbs.snr = np.full_like(sim_shape, np.nan)

        # Keep all frbs
        snr_mask = np.ones_like(frbs.snr, dtype=bool)

        # Have to loop over the observing times
        for i, t_min in enumerate(tqdm(times[:-1], desc='Pointings')):
            ra_p = survey.pointings[0][p_i % max_n_pointings]  # Loops around
            dec_p = survey.pointings[1][p_i % max_n_pointings]
            lst = lsts[i]
            t_max = times[i+1]

            # Which frbs are within the pointing time?
            time_mask = ((frbs.time >= t_min) & (frbs.time <= t_max))
            tm = np.any(time_mask, axis=1)  # Time mask (1D)

            # Of those, which are within the beam size?
            offset = go.separation(frbs.ra[tm], frbs.dec[tm], ra_p, dec_p)
            pm = (offset <= r)  # Position mask (1D)

            # Relevant FRBS
            tpm_1d = np.copy(tm)  # Time & position mask
            tnpm_1d = np.copy(tpm_1d)  # Time & not position mask
            tnpm_1d[tpm_1d] = ~pm
            tpm_1d[tpm_1d] = pm  # 1D
            tpm_2d = time_mask & tpm_1d[:, None]  # 2D
            tnpm_2d = time_mask & tnpm_1d[:, None]
            n_bursts = np.count_nonzero(tpm_2d, axis=1)  # Number of bursts

            # What's the intensity of them in the beam?
            int_pro, _ = survey.calc_int_pro(ra=frbs.ra[tpm_1d],
                                             dec=frbs.dec[tpm_1d],
                                             ra_p=ra_p,
                                             dec_p=dec_p,
                                             lst=lst,
                                             random_loc=False)

            # Ensure relevant masks
            mask = tpm_1d
            not_mask = tnpm_1d
            if frbs.s_peak.ndim == 2:
                mask = tpm_2d
                not_mask = tnpm_2d

            w_eff_mask = tpm_1d
            if frbs.w_eff.ndim == 2:
                w_eff_mask = tpm_2d

            # Apply intensities to those bursts' s_peak
            frbs.s_peak[mask] *= np.repeat(int_pro, n_bursts[n_bursts > 0])
            frbs.s_peak[not_mask] = np.nan

            # Ensure dimensionality is correct
            s_peak = frbs.s_peak[mask]
            w_eff = frbs.w_eff[w_eff_mask]
            if frbs.w_eff.ndim < frbs.s_peak.ndim:
                w_eff = np.repeat(frbs.w_eff, n_bursts)
            elif frbs.w_eff.ndim > frbs.s_peak.ndim:
                s_peak = np.repeat(frbs.s_peak, n_bursts)

            # Calculate fluence [Jy*ms]
            frbs.fluence[mask] = survey.calc_fluence(s_peak, w_eff)

            # Construct masks with correct dimension
            if frbs.w_arr.ndim == 2:
                w_arr = frbs.w_arr[tpm_2d]
            elif frbs.w_arr.ndim == 1:
                w_arr = frbs.w_arr[tpm_1d]

            # Ensure entries are repeated for the number of bursts
            s_peak = frbs.s_peak[mask]
            if frbs.w_arr.ndim < frbs.s_peak.ndim:
                w_arr = np.repeat(frbs.w_arr, n_bursts)
            elif frbs.w_arr.ndim > frbs.s_peak.ndim:
                s_peak = np.repeat(frbs.s_peak, n_bursts)

            # And system temperatures
            if frbs.T_sys.ndim == 1:
                T_sys = np.repeat(frbs.T_sys[tpm_1d], n_bursts[n_bursts > 0])
            elif frbs.T_sys.ndim == 0:
                T_sys = frbs.T_sys

            # Caculate Signal to Noise Ratio
            frbs.snr[mask] = survey.calc_snr(s_peak, w_arr, T_sys)

            # Add scintillation
            if self.scin:
                # Not been fully tested with repeaters
                # Might break due to differing dimensionality of dist and snr
                t_scat = frbs.t_scat[tpm_1d]
                dist_co = frbs.dist_co[tpm_1d]
                gl = frbs.gl[tpm_1d]
                gb = frbs.gb[tpm_1d]
                snr = frbs.snr[mask]
                new_snr = survey.calc_scint(t_scat, dist_co, gl, gb, snr)
                frbs.snr[mask] = np.repeat(new_snr, n_bursts[n_bursts > 0])

            # Don't keep those in time, but not in position
            snr_mask[not_mask] = False
            # Don't keep those in time, in position, but below S/N threshold
            snr_mask[mask] = (frbs.snr[mask] >= survey.snr_limit)

            if survey.strategy == 'regular':
                p_i += 1  # Head to next pointing

            # TODO: A 'follow-up' strategy could be implemented here
            # depending on the number of FRBs observed

        # TODO: Implement rates for sources vs bursts for repeater population

        # You should now have all frbs outside time
        # And frbs inside time, in position and above S/N
        frbs.apply(snr_mask)

        frbs.clean_up()

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
        self.rate.name = self.name
        self.rate.vol = self.rate.tot()
        self.rate.vol /= self.vol_co_max * (365.25*86400/self.time)

        # If oneoffs, you'll want to scale by the area and potentially time
        if not self.repeaters:
            self.rate.f_area = f_area
            self.rate.f_time = f_time

    def rates(self, scale_area=None, scale_time=None):
        """Scale the rates by the beam area and time."""
        r = self.rate

        if scale_area is None:
            scale_area = False
            if not self.repeaters:
                scale_area = True

        if scale_time is None:
            scale_time = False

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
