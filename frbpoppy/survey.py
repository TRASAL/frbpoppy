"""Class holding survey properties."""

from datetime import datetime, timedelta
from scipy.special import j1
import math
import numpy as np
import os
import pandas as pd

import frbpoppy.galacticops as go
from frbpoppy.log import pprint
from frbpoppy.paths import paths


class Survey:
    """
    Method containing survey parameters and functions.

    Args:
        name (str): Name of survey with which to observe population. Can
            either be a predefined survey present in frbpoppy or a path name to
            a new survey filename
        gain_pattern (str): Set gain pattern
        n_sidelobes (float): Number of sidelobes to include. Use 0.5 to
            cut beam at FWHM
        n_days (float): Time spent surveying [days]
        strategy (str): 'follow-up' or 'regular' (for RepeaterPopulation)

    """

    def __init__(self,
                 name,
                 gain_pattern='gaussian',
                 n_sidelobes=0.5,
                 n_days=1,
                 strategy='follow-up'):
        """Initializing."""
        # Set up parameters
        self.name = name
        self.gain_pattern = gain_pattern
        self.n_sidelobes = n_sidelobes
        self.n_days = n_days
        self.strategy = strategy
        self.beam_size = None
        self.pointings = None

        # Parse survey file
        self.read_survey_parameters()

        # Special treatment for perfect survey
        if self.name == 'perfect':
            self.gain_pattern = 'perfect'

        if self.transit is True:
            self.strategy = 'regular'

    def __str__(self):
        """Define how to print a survey object to a console."""
        s = 'Survey properties:'

        attributes = []
        for e in self.__dict__:
            attr = '\n\t{0:13.12}{1:.60}'.format(e, str(self.__dict__[e]))
            attributes.append(attr)

        s += ''.join(attributes)

        return s

    def read_survey_parameters(self):
        """Read in survey parameters."""
        # Read the survey file
        path = os.path.join(paths.surveys(), 'surveys.csv')
        df = pd.read_csv(path)
        df = df.set_index('survey')

        # Obtain the information relevant to a single survey
        survey = df[df.index == self.name].squeeze()

        # Parse parameters
        self.beta = survey['survey degradation factor']
        self.gain = survey['antenna gain (K/Jy)']
        self.t_obs = survey['integration time (s)']
        self.t_samp = survey['sampling time (ms)']
        self.T_rec = survey['receiver temperature (K)']
        self.central_freq = int(survey['centre frequency (MHz)'])
        self.bw = survey['bandwidth (MHz)']
        self.bw_chan = survey['channel bandwidth (MHz)']
        self.n_pol = survey['number of polarizations']
        self.beam_size_fwhm = survey['beam size (deg^2)']
        self.snr_limit = survey['signal-to-noise ratio [0-1]']
        self.max_w_eff = survey['maximum pulse width (ms)']
        self.latitude = survey['latitude (deg)']
        self.longitude = survey['longitude (deg)']
        self.transit = survey['transit telescope (bool)']
        self.ra_min = survey['minimum RA (deg)']
        self.ra_max = survey['maximum RA (deg)']
        self.dec_min = survey['minimum DEC (deg)']
        self.dec_max = survey['maximum DEC (deg)']
        self.gl_min = survey['minimum Galactic longitude (deg)']
        self.gl_max = survey['maximum Galactic longitude (deg)']
        self.gb_min = survey['minimum Galactic latitude (deg)']
        self.gb_max = survey['maximum Galactic latitude (deg)']
        self.up_time = survey['fractional uptime [0-1]']

    def in_region(self, frbs=None, ra=None, dec=None, gl=None, gb=None):
        """
        Check if the given frbs are within the survey region.

        Args:
            frbs (Frbs): Frbs of which to check whether in survey region
            ra, dec, gl, gb (float): Coordinates to use if not using frbs

        Returns:
            array: Boolean mask denoting whether frbs are within survey region

        """
        # Check whether to use user input or frbs coordinates
        if ra is None or dec is None:
            ra = frbs.ra
            dec = frbs.dec
        if gl is None or gb is None:
            gl = frbs.gl
            gb = frbs.gb

        # Create mask with False
        mask = np.ones_like(ra, dtype=bool)

        # Ensure in correct format
        gl[gl > 180.] -= 360.

        # Create region masks
        gl_limits = (gl > self.gl_max) | (gl < self.gl_min)
        gb_limits = (gb > self.gb_max) | (gb < self.gb_min)
        ra_limits = (ra > self.ra_max) | (ra < self.ra_min)
        dec_limits = (dec > self.dec_max) | (dec < self.dec_min)

        mask[gl_limits] = False
        mask[gb_limits] = False
        mask[ra_limits] = False
        mask[dec_limits] = False

        return mask

    def in_observation(self, frbs, strategy=None, t_stick=None,
                       pointings=None):
        """Whether an FRB has been detected using an observing strategy.

        Args:
            frbs (frbs): Frbs from a RepeaterPopulation
            strategy (str): Either 'follow-up' or 'regular'
            t_stick (float): Time in seconds to stick on same FRB if detected
            pointings (tuple): Tuple with arrays of ra, dec

        Returns:
            array: Mask with detections

        """
        if strategy is None:
            strategy = self.strategy
        if pointings is None:
            pointings = self.pointings

        # Create mask with False
        mask = np.zeros_like(frbs.time, dtype=bool)

        # Set up a list of pointings if not given
        if pointings is None:
            self.gen_pointings()

        # t_obs in fractional days
        t_obs = self.t_obs / 86400

        # Array with times of each pointing
        ts = np.arange(0, self.n_days+t_obs, t_obs)
        t_ranges = [(ts[i], ts[i+1]) for i in range(len(ts) - 1)]

        # Calculate size of detection region
        if self.beam_size is None:
            self.beam_size = self.beam_size_fwhm
        # Check whether the full sky
        if np.allclose(self.beam_size, 4*np.pi*(180/np.pi)**2):
            r = 180
        else:
            cos_r = (1 - (self.beam_size*np.pi)/(2*180**2))
            r = np.rad2deg(np.arccos(cos_r))

        i_p = 0  # Iterator for pointings
        if t_stick is None:
            t_stick = 1*t_obs  # Time to stick on pointing
        t_next_pointing = 0  # Time of next pointing
        n_p = len(self.pointings)  # Number of pointings

        for i_t, t_range in enumerate(t_ranges):
            # Split out time range
            t_min, t_max = t_range

            # Check whether FRBs are within pointing regions
            # The modulo means it will wrap back around to the first pointing
            ra, dec = self.pointings
            ra = ra[i_p % n_p]
            dec = dec[i_p % n_p]
            limit_pos = (go.separation(frbs.ra, frbs.dec, ra, dec) <= r)
            limit_pos = limit_pos[:, None]
            limit_time = ((frbs.time >= t_min) & (frbs.time <= t_max))
            limit = (limit_pos & limit_time)
            mask[limit] = True

            if strategy == 'follow-up':  # Follow up FRB
                if not np.any(limit):  # If no FRBs are detected
                    if t_min >= t_next_pointing:  # If not sticking to field
                        i_p += 1  # Go to next pointing
                else:  # If there are FRBs
                    # Set time after which you can head to the next pointing
                    t_next_pointing = t_min + t_stick
            elif strategy == 'regular':  # Just stick to your regular schedule
                i_p += 1

        return mask

    def gen_pointings(self):
        """Generate pointings.

        Returns:
            ra, dec: Arrays with ra, dec of pointings

        """
        n_p = int(self.n_days*24*60*60 / self.t_obs)  # Number of pointings
        if self.transit:
            self.pointings = self.gen_transit_pointings(n_p,
                                                        self.latitude,
                                                        self.longitude,
                                                        self.t_obs)
        else:
            self.pointings = self.gen_tracking_pointings(n_p)

    def gen_transit_pointings(self, n_gen, lat=None, lon=None, t_obs=None):
        """Generate RA and Dec pointing coordinates for a transit telescope.

        Args:
            n_gen (int): Number of pointings wanted.
            lat (float): Latitude of telescope.
            lon (float): Longitude of telescope (minus for west).
            t_obs (type): Time for each pointing.

        Returns:
            ra, dec: Numpy arrays with coordinates of pointings

        """
        if lat is None:
            lat = self.latitude
            lon = self.longitude
            t_obs = self.t_obs

        # Pointings are randomly placed between the year 2000 and 2100
        date_min = go.random_date(datetime(2000, 1, 1), datetime(2100, 1, 1))
        date_max = date_min + timedelta(seconds=int(t_obs*n_gen))
        time_delta = np.timedelta64(t_obs, 's')
        times = np.arange(date_min, date_max, time_delta, dtype='datetime64')

        ra = go.datetime_to_gmst(times) + lon
        dec = np.ones(n_gen)*lat

        return ra, dec

    def gen_tracking_pointings(self, n_gen):
        """Generate RA, Dec pointing coordinates for a tracking telescope.

        Uses a try-accept algorithm to generate pointings with a survey. For
        details on the sunflower algoritm used to distribute pointings see
        See https://stackoverflow.com/a/44164075/11922471

        Args:
            n_gen (int): Number of pointings.

        Returns:
            ra, dec: Numpy arrays with coordinates of pointings

        """
        def sample(n=1000, random=True):
            """Length of sunflower-like chain to sample."""
            indices = np.arange(0, n, dtype=float) + 0.5

            # Generate coordinates using golden ratio
            ra = (np.pi * (1 + 5**0.5) * indices) % (2*np.pi)
            dec = np.arccos(1 - 2*indices/n) - 0.5*np.pi

            if random:
                # Start from a random point rather than the south pole
                phi = np.random.uniform(-np.pi, np.pi, 1)  # Random ra
                theta = np.random.uniform(-np.pi/2, np.pi/2, 1)  # Random dec

                # Shift ra from -180 to 180
                ra[ra > np.pi] -= 2*np.pi

                # Save on line length
                sin = np.sin
                cos = np.cos
                arcsin = np.arcsin

                y = sin(ra)
                x = np.tan(dec)*sin(theta) + cos(ra)*cos(theta)
                lon = np.arctan2(y, x) - phi
                lat = arcsin(cos(theta)*sin(dec)-cos(ra)*sin(theta)*cos(dec))

                # Shift ra back to 0-360
                ra[ra < 0] += 2*np.pi
                ra = lon % (2*np.pi)
                dec = lat

            # To degrees
            ra = np.rad2deg(ra)
            dec = np.rad2deg(dec)

            return ra, dec

        def accept(ra, dec):
            gl, gb = go.radec_to_lb(ra, dec, frac=True)
            return self.in_region(ra=ra, dec=dec, gl=gl, gb=gb)

        # Accept-Reject
        n = n_gen
        ra, dec = sample(n)
        mask = accept(ra, dec)

        # While there are not enough valid points, keep generating
        while sum(mask) < n_gen:
            n += (n_gen - sum(mask))
            ra, dec = sample(n)
            mask = accept(ra, dec)
        else:
            # Only select valid points
            ra = ra[mask]
            dec = dec[mask]

            # Evenly sample arrays
            idx = np.round(np.linspace(0, len(ra) - 1, n_gen)).astype(int)
            ra = ra[idx]
            dec = dec[idx]

        return ra, dec

    def max_offset(self, x):
        """Calculate the maximum offset of an FRB in an Airy disk.

        Args:
            x (int): Maximum sidelobe wanted

        """
        # Null points of kasin for allow a number of sidelobes
        kasin_nulls = [3.831706, 7.015587, 10.173468, 13.323692, 16.47063,
                       19.615859, 22.760084, 25.903672, 29.046829, 32.18968,
                       35.332308, 38.474766]

        # Allow for cut at FWHM
        if x == 0.5:
            return 1

        try:
            arcsin = math.asin(self.fwhm*kasin_nulls[x]/(60*180))
        except ValueError:
            m = f'Beamsize including sidelobes would be larger than sky \n'
            A = (90/kasin_nulls[x])**2*math.pi
            m += f'Ensure beamsize is smaller than {A}'
            raise ValueError(m)

        return 2/self.fwhm * 60*180/math.pi * arcsin

    def intensity_profile(self, shape=1, dimensions=2, rep_loc='random'):
        """Calculate intensity profile.

        Args:
            shape (tuple): Usually the shape of frbs.s_peak
            dimensions (int): Use a 2D beampattern or a 1D one.
            rep_loc (str): 'same' or 'random'. Whether repeaters are observed
                in the same spot or a different spot

        Returns:
            array, array: intensity profile, offset from beam [arcmin]

        """
        # Calculate Full Width Half Maximum from beamsize
        cos_r = (1 - (self.beam_size_fwhm*np.pi)/(2*180**2))
        offset = np.rad2deg(np.arccos(cos_r)) * 60  # [arcmin]
        self.fwhm = 2*offset  # [arcmin]

        if rep_loc == 'same':
            r = np.random.random(shape[0])
        else:
            r = np.random.random(shape)

        if dimensions == 2:  # 2D
            offset *= np.sqrt(r)
        elif dimensions == 1:  # 1D
            offset *= r

        # Allow for a perfect beam pattern in which all is detected
        if self.gain_pattern == 'perfect':
            int_pro = np.ones(shape)
            self.beam_size = self.beam_size_fwhm
            return int_pro, offset

        # Formula's based on 'Interferometry and Synthesis in Radio
        # Astronomy' by A. Richard Thompson, James. M. Moran and
        # George W. Swenson, JR. (Second edition), around p. 15

        max_offset = self.max_offset(self.n_sidelobes)
        self.beam_size = math.pi*(self.fwhm/2*max_offset/60)**2  # [sq degrees]

        if self.gain_pattern == 'gaussian':
            # Set the maximum offset equal to the null after a sidelobe
            # I realise this pattern isn't an airy, but you have to cut
            # somewhere
            offset *= max_offset
            alpha = 2*math.sqrt(math.log(2))
            int_pro = np.exp(-(alpha*offset/self.fwhm)**2)
            return int_pro, offset

        elif self.gain_pattern == 'airy':
            # Set the maximum offset equal to the null after a sidelobe
            offset *= max_offset
            c = 299792458
            conv = math.pi/(60*180)  # Conversion arcmins -> radians
            eff_diam = c/(self.central_freq*1e6*conv*self.fwhm)
            a = eff_diam/2  # Effective radius of telescope
            lamda = c/(self.central_freq*1e6)
            ka = (2*math.pi*a/lamda)
            kasin = ka*np.sin(offset*conv)
            int_pro = 4*(j1(kasin)/kasin)**2
            return int_pro, offset

        elif self.gain_pattern in ['parkes', 'apertif']:

            place = paths.models() + f'/beams/{self.gain_pattern}.npy'
            beam_array = np.load(place)
            b_shape = beam_array.shape
            ran_x = np.random.randint(0, b_shape[0], shape)
            ran_y = np.random.randint(0, b_shape[1], shape)
            int_pro = beam_array[ran_x, ran_y]
            offset = np.sqrt((ran_x-b_shape[0]/2)**2 + (ran_y-b_shape[1]/2)**2)

            # Scaling factors to correct for pixel scale
            if self.gain_pattern == 'apertif':  # 1 pixel = 0.94'
                offset *= 240/256  # [arcmin]
                self.beam_size = 25.
            if self.gain_pattern == 'parkes':  # 1 pixel = 54"
                offset *= 0.9  # [arcmin]
                self.beam_size = 9.

            return int_pro, offset

        else:
            pprint(f'Gain pattern "{self.gain_pattern}" not recognised')

    def dm_smear(self, frbs):
        """
        Calculate delay in pulse across a channel due to dm smearing.

        Formula's based on 'Handbook of Pulsar Astronomy" by Duncan Lorimer
        & Michael Kramer, section A2.4. Note the power of the forefactor has
        changed due to the central frequency being given in MHz.

        Args:
            frbs (FRBs): FRB object with a dm attribute

        Returns:
            t_dm (array): Time of delay [ms] at central band frequency

        """
        return 8.297616e6 * self.bw_chan * frbs.dm * (self.central_freq)**-3

    def calc_scat(self, dm):
        """Calculate scattering timescale for FRBs.

        Offset according to Lorimer et al. (doi:10.1093/mnrasl/slt098)

        Args:
            dm (array): Dispersion Measure

        Returns:
            array: Scattering timescales [ms]

        """
        freq = self.central_freq
        return go.scatter_bhat(dm, scindex=-3.86, offset=-9.5, freq=freq)

    def calc_Ts(self, frbs):
        """Set temperatures for frbs."""
        # Special treatment for a perfect telescope
        if self.name.startswith('perfect'):
            T_sky = 0
            T_sys = self.T_rec
        else:
            T_sky = self.calc_T_sky(frbs)
            T_sys = self.T_rec + T_sky

        return T_sky, T_sys

    def calc_T_sky(self, frbs):
        """
        Calculate the sky temperature from the Haslam table.

        Afterwards scale to the survey frequency. The temperature sky map is
        given in the weird units of HealPix and despite looking up info on this
        coordinate system, I don't have the foggiest idea of how to transform
        these to galactic coordinates. I have therefore directly copied the
        following code from psrpoppy in the assumption Sam Bates managed to
        figure it out.

        Args:
            frbs (FRBs): Needed for coordinates
        Returns:
            array: Sky temperature [K]
        """
        T_sky_list = go.load_T_sky()

        # ensure l is in range 0 -> 360
        B = frbs.gb
        L = np.copy(frbs.gl)
        L[L < 0.] += 360

        # convert from l and b to list indices
        j = B + 90.5
        j[j > 179] = 179

        nl = L - 0.5
        nl[L < 0.5] = 359
        i = nl / 4.

        index = 180*i.astype(int) + j.astype(int)
        T_sky_haslam = np.take(T_sky_list, index).astype(np.float32)

        # scale temperature
        # Assuming dominated by syncrotron radiation
        T_sky = T_sky_haslam * (self.central_freq/408.0)**(-2.6)

        return T_sky

    def calc_s_peak(self, frbs, f_low=100e6, f_high=10e9):
        """
        Calculate the mean spectral flux density.

        Following Lorimer et al, 2013, eq. 9., at the central frequency
        of the survey.

        Args:
            frbs (FRBs): FRBs
            f_low (float): Source emission lower frequency limit [Hz].
            f_high (float): Source emission higher frequency limit [Hz].

        Returns:
            array: Mean spectral flux density [Jy]

        """
        # Limits observing bandwidth (as seen in rest frame source)
        f_1 = (self.central_freq - 0.5*self.bw)
        f_1 *= 1e6  # MHz -> Hz
        f_2 = (self.central_freq + 0.5*self.bw)
        f_2 *= 1e6  # MHz -> Hz

        # Spectral index
        si = frbs.si
        lum_bol = frbs.lum_bol
        if lum_bol.ndim == si.ndim:
            pass
        elif lum_bol.ndim > si.ndim:
            si = si[:, None]
        elif lum_bol.ndim < si.ndim:
            lum_bol = lum_bol[:, None]

        sp = si + 1
        sm = si - 1

        freq_frac = (f_2**sp - f_1**sp) / (f_2 - f_1)

        z = frbs.z

        # Convert distance in Gpc to 10^25 metres
        dist = frbs.dist_co * 3.085678
        dist = dist.astype(np.float64)

        # Convert luminosity in 10^7 Watts such that s_peak will be in Janskys
        lum = lum_bol * 1e-31
        lum = lum.astype(np.float64)

        if lum.ndim > 1:
            z = z[:, None]
            dist = dist[:, None]

        nom = lum * (1+z)**sm * freq_frac
        den = 4*np.pi*dist**2 * (f_high**sp - f_low**sp)

        if nom.ndim == den.ndim:
            pass
        elif nom.ndim > 1:
            den = den[:, None]

        s_peak = nom/den

        # Add degradation factor due to pulse broadening (see Connor 2019)
        w_frac = (frbs.w_arr / frbs.w_eff)

        # Dimension conversions
        if w_frac.ndim == s_peak.ndim:
            pass
        elif s_peak.ndim == 1:
            s_peak = s_peak[:, None]
        elif w_frac.ndim == 1:
            w_frac = w_frac[:, None]

        s_peak = s_peak * w_frac

        return s_peak.astype(np.float32)

    def calc_w_eff(self, frbs):
        """Calculate effective pulse width [ms].

        From Narayan (1987, DOI: 10.1086/165442), and also Cordes & McLaughlin
        (2003, DOI: 10.1086/378231). For details see p. 30 of Emily Petroff's
        thesis (2016), found here: http://hdl.handle.net/1959.3/417307

        Args:
            frbs (FRBs): FRBs for which to calculate effective pulse width.

        Returns:
            array: Effective pulse width [ms]

        """
        w_arr = frbs.w_arr
        t_dm = frbs.t_dm
        t_scat = frbs.t_scat

        if frbs.w_arr.ndim > 1:
            t_dm = t_dm[:, None]
            if isinstance(frbs.t_scat, np.ndarray):
                t_scat = t_scat[:, None]

        return np.sqrt(w_arr**2 + t_dm**2 + t_scat**2 + self.t_samp**2)

    def calc_snr(self, frbs):
        """
        Caculate the SNR of several frbs.

        Radiometer equation for single pulse (Dewey et al., 1984), but
        adapted to allow for a degradation factor reducing the peak flux
        as a pulse is stretched

        Args:
            frbs (FRBs): FRBs of which to calculate the signal to noise

        Returns:
            array: Signal to noise ratio based on the radiometer
                equation for a single pulse.

        """
        sp = frbs.s_peak
        w_arr = frbs.w_arr
        T_sys = frbs.T_sys

        # Dimension check
        if sp.ndim == w_arr.ndim:
            pass
        elif sp.ndim == 1:
            sp = sp[:, None]
        elif w_arr.ndim == 1:
            w_arr = w_arr[:, None]

        if (sp.ndim or w_arr.ndim) > 1:
            if isinstance(T_sys, np.ndarray):
                T_sys = T_sys[:, None]

        snr = sp*self.gain*np.sqrt(self.n_pol*self.bw*w_arr*1e3)
        snr /= (T_sys * self.beta)

        return snr

    def calc_fluence(self, frbs):
        """Calculate fluence [Jy*ms]."""
        if frbs.w_eff.ndim == frbs.s_peak.ndim:
            return frbs.s_peak * frbs.w_eff
        elif frbs.w_eff.ndim == 1:
            return frbs.s_peak * frbs.w_eff[:, None]
        elif frbs.s_peak.ndim == 1:
            return frbs.s_peak[:, None] * frbs.w_eff

    def calc_scint(self, frbs):
        """
        Calculate scintillation effect on the signal to noise ratio.

        (Rather than adapting the flux, as the snr can change per survey
        attempt). Formulas based on 'Handbook of Pulsar Astronomy" by Duncan
        Lorimer & Michael Kramer, section 4.2. Test this before applying - no
        rigorous testing has been applied to this.

        Args:
            frbs (FRBs): FRBs

        Returns:
            array: Signal to noise ratio modulation factors for scintillation

        """
        # Calculate scattering
        if type(frbs.t_scat) is not np.ndarray:
            frbs.t_scat = self.calc_scat(frbs.dm)

        # Convert to seconds
        frbs.t_scat /= 1000.

        # Decorrelation bandwidth (eq. 4.39)
        decorr_bw = 1.16/(2*math.pi*frbs.t_scat)
        # Convert to MHz
        decorr_bw /= 1e6

        # Scintillation strength (eq. 4.33)
        u = np.sqrt(self.central_freq / decorr_bw)

        m = np.zeros_like(u)

        # Strong scintillation
        strong = (u < 1)
        m[strong] = np.sqrt(u[strong]**(5/3))   # (eq. 4.35)

        # Weak scintillation

        # Refractive scintillation (eq. 4.47)
        m_riss = u**-(1/3)

        # Taking the average kappa value
        kappa = 0.15

        t_diss, decorr_bw = go.ne2001_scint_time_bw(frbs.dist_co,
                                                    frbs.gl,
                                                    frbs.gb,
                                                    self.central_freq)

        # Following Cordes and Lazio (1991) (eq. 4.43)
        n_t = np.ones_like(t_diss)
        n_t[~np.isnan(t_diss)] = 1 + kappa * self.t_obs / t_diss

        n_f = np.ones_like(decorr_bw)
        n_f[~np.isnan(decorr_bw)] = 1 + kappa * self.bw / decorr_bw

        # Diffractive scintillation (eq. 4.41)
        m_diss = 1 / np.sqrt(n_t * n_f)

        # (eq. 4.48)
        weak = (u >= 1)
        m[weak] = np.sqrt(m_diss**2 + m_riss**2 + m_diss*m_riss)

        # Distribute the scintillation according to gaussian distribution
        snr = np.random.normal(frbs.snr, m*frbs.snr)

        return snr

    def calc_fluence_limit(self, w_eff=None):
        """Calculate the fluence limit.

        Read Keane, Petroff (2015) for more details on how this is calculated.

        Args:
            w_eff: Pulse width at which to calculate the fluence limit [ms].
                Default sets this to be at the maximum searched pulse width.

        Returns:
            float: Fluence limit for a maximum pulse width burst [Jy ms]

        """
        if not w_eff:
            w_eff = self.max_w_eff

        # Line of constant S/N
        self.s_peak_limit = self.snr_limit*self.T_rec*self.beta
        self.s_peak_limit /= self.gain*math.sqrt(self.n_pol*self.bw*1e3)

        # Line of constant fluence
        self.fluence_limit = self.s_peak_limit / math.sqrt(w_eff)
        self.fluence_limit *= w_eff

        return self.fluence_limit
