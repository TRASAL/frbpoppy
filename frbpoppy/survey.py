"""Class holding survey properties."""

from scipy.special import j1
import numpy as np
import os
import pandas as pd

import frbpoppy.galacticops as go
import frbpoppy.pointings as pointings
from frbpoppy.paths import paths


class Survey:
    """
    Method containing survey parameters and functions.

    Args:
        name (str): Name of survey with which to observe population. Can
            either be a predefined survey present in frbpoppy or a path name to
            a new survey filename
        n_days (float): Time spent surveying [days]

    """

    def __init__(self,
                 name='perfect',
                 n_days=1):
        """Initializing."""
        # Set up parameters
        self.name = name
        self.n_days = n_days

        # Calculate some at a later stage
        self.beam_size = None
        self.beam_array = None
        self.pointings = None

        # Parse survey file
        self.read_survey_parameters()

        # Special treatment for perfect survey
        if self.name.startswith('perfect'):
            self.beam_pattern = 'perfect'

        # Un-used
        self.strategy = 'regular'  # 'follow-up' not implemented
        if self.mount_type is 'transit':
            # Transit telescopes can't follow-up
            self.strategy = 'regular'

        # Set beam properties
        self.set_beam()

        # Set pointing properties
        self.set_pointings(mount_type=self.mount_type)

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
        self.mount_type = survey['mount type']
        self.ra_min = survey['minimum RA (deg)']
        self.ra_max = survey['maximum RA (deg)']
        self.dec_min = survey['minimum DEC (deg)']
        self.dec_max = survey['maximum DEC (deg)']
        self.gl_min = survey['minimum Galactic longitude (deg)']
        self.gl_max = survey['maximum Galactic longitude (deg)']
        self.gb_min = survey['minimum Galactic latitude (deg)']
        self.gb_max = survey['maximum Galactic latitude (deg)']

    def calc_sky_radius(self, area=None):
        """Determine the radius along the sky of an area in sq. degrees."""
        # Calculate size of detection region
        if area is None:
            area = self.beam_size_fwhm

        # Check whether the full sky
        if np.allclose(area, 4*np.pi*(180/np.pi)**2):
            r = 180
        else:
            cos_r = (1 - (area*np.pi)/(2*180**2))
            r = np.rad2deg(np.arccos(cos_r))

        return r

    def in_region(self, ra, dec, gl, gb):
        """
        Check if the given frbs are within the survey region.

        Args:
            ra, dec, gl, gb (float): Coordinates to check whether in region

        Returns:
            array: Boolean mask denoting whether frbs are within survey region

        """
        # Create mask with False if not in region
        mask = go.in_region(ra, dec, gl, gb,
                            ra_min=self.ra_min, ra_max=self.ra_max,
                            dec_min=self.dec_min, dec_max=self.dec_max,
                            gl_min=self.gl_min, gl_max=self.gl_max,
                            gb_min=self.gb_min, gb_max=self.gb_max)

        return mask

    def set_pointings(self, mount_type='tracking', n_pointings=None, ra=None,
                      dec=None):
        """Set pointing properties."""
        # If you want to use your own pointings
        if ra and dec:
            self.n_pointings = len(ra)
            self.point_func = lambda: (ra, dec)
            return

        # Calculate the number of pointings needed
        self.n_pointings = n_pointings
        if n_pointings is None:
            self.n_pointings = int(self.n_days*24*60*60 / self.t_obs)

        if mount_type == 'transit':
            self.point_func = lambda: pointings.transit(self.n_pointings,
                                                        self.latitude,
                                                        self.longitude,
                                                        self.t_obs)
        else:
            self.point_func = lambda: pointings.tracking(self.n_pointings,
                                                         ra_min=self.ra_min,
                                                         ra_max=self.ra_max,
                                                         dec_min=self.dec_min,
                                                         dec_max=self.dec_max,
                                                         gl_min=self.gl_min,
                                                         gl_max=self.gl_max,
                                                         gb_min=self.gb_min,
                                                         gb_max=self.gb_max)

    def gen_pointings(self):
        """Generate pointings."""
        self.pointings = self.point_func()

    def max_offset(self, n):
        """Calculate the maximum offset of an FRB in an Airy disk.

        Args:
            n (int): Maximum number of wanted sidelobes

        """
        # Null points of kasin for allow a number of sidelobes
        kasin_nulls = [3.831706, 7.015587, 10.173468, 13.323692, 16.47063,
                       19.615859, 22.760084, 25.903672, 29.046829, 32.18968,
                       35.332308, 38.474766]

        # Allow for cut at FWHM
        if n == 0.5:
            return 1

        arcsin = np.arcsin(self.fwhm*kasin_nulls[n]/180)
        if np.isnan(arcsin):
            m = f'Beamsize including sidelobes would be larger than sky \n'
            A = (90/kasin_nulls[n])**2*np.pi
            m += f'Ensure beamsize is smaller than {A}'
            raise ValueError(m)

        return 2/self.fwhm * 180/np.pi * arcsin

    def calc_random_int_pro(self, shape=1):
        """Calculate the intensity profile in random spots of a beam.

        Args:
            shape (tuple): Usually the shape of frbs.s_peak

        Returns:
            array, array: intensity profile, offset from beam [degree]

        """
        # Calculate Full Width Half Maximum from beamsize
        offset = self.calc_sky_radius()
        self.fwhm = 2*offset  # The diameter [deg]

        # Take a random location in the 2D beampattern
        offset *= np.sqrt(np.random.random(shape))

        # Allow for a perfect beam pattern in which all is detected
        if self.beam_pattern.startswith('perfect'):
            int_pro = np.ones(shape)
            if self.beam_size is None:
                self.beam_size = self.beam_size_fwhm
            return int_pro, offset

        # Formula's based on 'Interferometry and Synthesis in Radio
        # Astronomy' by A. Richard Thompson, James. M. Moran and
        # George W. Swenson, JR. (Second edition), around p. 15

        max_offset = self.max_offset(self.n_sidelobes)

        if self.beam_size is None:
            self.beam_size = np.pi*(self.fwhm/2*max_offset)**2  # [sq deg]

        if self.beam_pattern == 'gaussian':
            # Set the maximum offset equal to the null after a sidelobe
            # I realise this pattern isn't an airy, but you have to cut
            # somewhere
            offset *= max_offset
            alpha = 2*np.sqrt(np.log(2))
            int_pro = np.exp(-(alpha*offset/self.fwhm)**2)
            return int_pro, offset

        elif self.beam_pattern == 'airy':
            # Set the maximum offset equal to the null after a sidelobe
            offset *= max_offset
            c = 299792458
            conv = np.pi/180  # Conversion degrees -> radians
            eff_diam = c/(self.central_freq*1e6*conv*self.fwhm)
            a = eff_diam/2  # Effective radius of telescope
            lamda = c/(self.central_freq*1e6)
            ka = (2*np.pi*a/lamda)
            kasin = ka*np.sin(offset*conv)
            int_pro = 4*(j1(kasin)/kasin)**2
            return int_pro, offset

        elif self.beam_array is not None:
            b_shape = self.beam_array.shape
            ran_x = np.random.randint(0, b_shape[0], shape)
            ran_y = np.random.randint(0, b_shape[1], shape)
            int_pro = self.beam_array[ran_x, ran_y]
            offset = np.sqrt((ran_x-b_shape[0]/2)**2 + (ran_y-b_shape[1]/2)**2)
            offset *= self.pixel_scale  # Correct for pixel scale
            return int_pro, offset

        else:
            m = f'Beam pattern "{self.beam_pattern}" not recognised'
            raise ValueError(m)

    def calc_fixed_int_pro(self, ra, dec, ra_p, dec_p, lst, test=False):
        """Calculate intensity profile for fixed location in beam.

        Args:
            ra (array): Right ascension of objects [deg]
            dec (array): Declination of objects [deg]
            ra_p (float): Right ascension of pointing [deg]
            dec_p (float): Declination of pointing [deg]
            lst (float): Local Sidereal Time [deg]
            test (bool): For testing

        Returns:
            type: Description of returned object.

        """
        # Weed out perfect beam
        if self.beam_pattern.startswith('perfect'):
            if self.beam_size is None:
                self.beam_size = self.beam_size_fwhm
            return np.ones(len(ra))

        # Convert input decimal degrees to radians
        ra = np.deg2rad(ra)
        dec = np.deg2rad(dec)
        args = [ra_p, dec_p, lst, self.latitude, self.pixel_scale]

        for a in args:
            if a is None:
                raise ValueError('Missing required input')

        ra_p, dec_p, lst, lat, pixel_scale = np.deg2rad(args)

        if self.mount_type == 'equatorial':
            # Convert input coordinates to offset in ra and dec
            dx, dy = go.coord_to_offset(ra_p, dec_p, ra, dec)
        elif self.mount_type in ('azimuthal', 'transit'):
            # Convert input right ascension to hour angle
            ha = lst - ra
            ha_p = lst - ra_p
            # Convert ha, dec to az, alt
            az, alt = go.hadec_to_azalt(ha, dec, lat)
            az_p, alt_p = go.hadec_to_azalt(ha_p, dec_p, lat)
            # Convert to offset
            dx, dy = go.coord_to_offset(az_p, alt_p, az, alt)
        else:
            raise ValueError(f'Invalid mount type: {self.mount_type}')

        # Convert offsets dx, dy to pixel in beam pattern (round)
        dx_px = (np.round(dx / pixel_scale)).astype(int)
        dy_px = (np.round(dy / pixel_scale)).astype(int)
        ny, nx = self.beam_array.shape
        x = (nx/2 + dx_px).astype(int)
        y = (ny/2 + dy_px).astype(int)

        # Get the value at this pixel (zero if outside beam pattern)
        m = self.beam_array.shape[0]
        outside = ((x <= 0) | (x >= m) | (y <= 0) | (y >= m))
        x[outside] = 0  # Nans don't work in int arrays
        y[outside] = 0

        intensity = self.beam_array[y, x]
        intensity[(x == 0) & (y == 0)] = 0

        return intensity

    def set_beam(self, model='perfect', random_loc=True, n_sidelobes=0.5):
        """Set intensity profile.

        Set properties for int pro

        Args:
            model (str): Beam pattern. Choice from 'apertif', 'parkes',
                'chime', 'gaussian', 'airy'.
            random_loc (bool): Whether to calculate the location of each burst
                or place them randomly in the beam.
        if model in ('gaussian', 'airy'):
            n_sidelobes (int): Number of sidelobes to include.

        """
        # Set up beam properties
        self.beam_pattern = model
        self.n_sidelobes = n_sidelobes

        # Set up predicted beam models
        if model in ('apertif', 'parkes', 'chime', 'gaussian', 'airy'):
            place = paths.models() + f'/beams/{self.beam_pattern}.npy'
            self.beam_array = np.load(place)

            if model == 'apertif':
                self.pixel_scale = 0.94/60  # Degrees per pixel [deg]
                self.beam_size = 25.  # [sq deg]
            elif model == 'parkes':
                self.pixel_scale = 54/3600  # Degrees per pixel [deg]
                self.beam_size = 9.  # [sq deg]
            elif model == 'chime':
                self.pixel_scale = 0.08  # Degrees per pixel [deg]
                self.beam_size = 80*80  # [sq deg]
            elif model == 'gaussian':
                r = self.calc_sky_radius()
                self.pixel_scale = r / 95  # Degrees per pixel [deg]
                self.beam_size = (self.pixel_scale*self.beam_array.shape[0])**2
            elif model == 'airy':
                r = self.calc_sky_radius()
                self.pixel_scale = r / 31  # Degrees per pixel [deg]
                self.beam_size = (self.pixel_scale*self.beam_array.shape[0])**2
        elif not model.startswith('perfect'):
            raise ValueError('Beam model input not recognised.')

        if random_loc or model.startswith('perfect'):
            self.beam_func = self.calc_random_int_pro
        else:
            self.beam_func = self.calc_fixed_int_pro

    def calc_dm_smear(self, dm):
        """
        Calculate delay in pulse across a channel due to dm smearing.

        Formula's based on 'Handbook of Pulsar Astronomy" by Duncan Lorimer
        & Michael Kramer, section A2.4. Note the power of the forefactor has
        changed due to the central frequency being given in MHz.

        Args:
            dm (array): Dispersion measure [pc/cm^3]

        Returns:
            t_dm (array): Time of delay [ms] at central band frequency

        """
        return 8.297616e6 * self.bw_chan * dm * (self.central_freq)**-3

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

    def calc_Ts(self, gl, gb):
        """Set temperatures for frbs.

        Args:
            gl (array): Galactic longitude [deg]
            gb (array): Galactic latitude [deg]

        Returns:
            T_sky, T_sys [K]

        """
        # Special treatment for a perfect telescope
        if self.name.startswith('perfect'):
            T_sky = 0
            T_sys = self.T_rec
        else:
            T_sky = self.calc_T_sky(gl, gb)
            T_sys = self.T_rec + T_sky

        return T_sky, T_sys

    def calc_T_sky(self, gl, gb):
        """
        Calculate the sky temperature from the Haslam table.

        Afterwards scale to the survey frequency. The temperature sky map is
        given in the weird units of HealPix and despite looking up info on this
        coordinate system, I don't have the foggiest idea of how to transform
        these to galactic coordinates. I have therefore directly copied the
        following code from psrpoppy in the assumption Sam Bates managed to
        figure it out.

        Args:
            gl (array): Galactic longitude [deg]
            gb (array): Galactic latitude [deg]
        Returns:
            array: Sky temperature [K]

        """
        T_sky_list = go.load_T_sky()

        # ensure l is in range 0 -> 360
        B = gb
        L = np.copy(gl)
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

    def calc_s_peak(self, si, lum_bol, z, dist_co, w_arr, w_eff,
                    f_low=100e6, f_high=10e9):
        """
        Calculate the mean spectral flux density.

        Following Lorimer et al, 2013, eq. 9., at the central frequency
        of the survey.

        Args:
            si (array): Spectral index
            lum_bol (array): Bolometric luminosity within emission band
            z (array): Redshift
            dist_co (array): Comoving distance [Gpc]
            w_arr (array): Pulse width at Earth [ms]
            w_eff (array): Pulse width at point of detection [ms]
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
        if lum_bol.ndim == si.ndim:
            pass
        elif lum_bol.ndim > si.ndim:
            si = si[:, None]
        elif lum_bol.ndim < si.ndim:
            lum_bol = lum_bol[:, None]

        sp = si + 1
        sm = si - 1

        freq_frac = (f_2**sp - f_1**sp) / (f_2 - f_1)

        # Convert distance in Gpc to 10^25 metres
        dist = dist_co * 3.085678
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
        w_frac = (w_arr / w_eff)

        # Dimension conversions
        if w_frac.ndim == s_peak.ndim:
            pass
        elif s_peak.ndim == 1:
            s_peak = s_peak[:, None]
        elif w_frac.ndim == 1:
            w_frac = w_frac[:, None]

        s_peak = s_peak * w_frac

        return s_peak.astype(np.float32)

    def calc_w_eff(self, w_arr, t_dm, t_scat):
        """Calculate effective pulse width [ms].

        From Narayan (1987, DOI: 10.1086/165442), and also Cordes & McLaughlin
        (2003, DOI: 10.1086/378231). For details see p. 30 of Emily Petroff's
        thesis (2016), found here: http://hdl.handle.net/1959.3/417307

        Args:
            w_arr (array): Pulse width [ms]
            t_dm (array): Time delay due to dispersion [ms]
            t_scat (array): Time delay due to scattering [ms]

        Returns:
            array: Effective pulse width [ms]

        """
        if w_arr.ndim > 1:
            t_dm = t_dm[:, None]
            if isinstance(t_scat, np.ndarray):
                t_scat = t_scat[:, None]

        return np.sqrt(w_arr**2 + t_dm**2 + t_scat**2 + self.t_samp**2)

    def calc_snr(self, s_peak, w_arr, T_sys):
        """
        Caculate the SNR of several frbs.

        Radiometer equation for single pulse (Dewey et al., 1984), but
        adapted to allow for a degradation factor reducing the peak flux
        as a pulse is stretched

        Args:
            s_peak (array): Peak flux [Jy]
            w_arr (array): Pulse width at Earth [ms]
            T_sys (array): System temperature [K]

        Returns:
            array: Signal to noise ratio based on the radiometer
                equation for a single pulse.

        """
        # Dimension check
        if s_peak.ndim == w_arr.ndim:
            pass
        elif s_peak.ndim == 1:
            s_peak = s_peak[:, None]
        elif w_arr.ndim == 1:
            w_arr = w_arr[:, None]

        if (s_peak.ndim or w_arr.ndim) > 1:
            if isinstance(T_sys, np.ndarray):
                T_sys = T_sys[:, None]

        snr = s_peak*self.gain*np.sqrt(self.n_pol*self.bw*w_arr*1e3)
        snr /= (T_sys * self.beta)

        return snr

    def calc_fluence(self, s_peak, w_eff):
        """Calculate fluence [Jy ms].

        Args:
            s_peak (array): Peak flux
            w_eff (array): Effective pulse width

        Returns:
            array: Fluence in Jy*ms

        """
        if w_eff.ndim == s_peak.ndim:
            return s_peak * w_eff
        elif w_eff.ndim == 1:
            return s_peak * w_eff[:, None]
        elif s_peak.ndim == 1:
            return s_peak[:, None] * w_eff

    def calc_scint(self, t_scat, dist_co, gl, gb, snr):
        """
        Calculate scintillation effect on the signal to noise ratio.

        (Rather than adapting the flux, as the snr can change per survey
        attempt). Formulas based on 'Handbook of Pulsar Astronomy" by Duncan
        Lorimer & Michael Kramer, section 4.2. Test this before applying - no
        rigorous testing has been applied to this.

        Args:
            t_scat (array): Scattering timescale for frbs [ms]
                (can use self.calc_scat to calculate this)
            dist_co (array): Comoving distance array [Gpc]
            gl (array): Galactic longitude [deg]
            gb (array): Galactic latitude [deg]
            snr (array): Signal to Noise array to modify

        Returns:
            array: Signal to noise ratio modulation factors for scintillation

        """
        if not isinstance(t_scat, np.ndarray):
            m = 'Please ensure t_scat is calculated using self.calc_scat'
            raise ValueError(m)

        # Convert to seconds
        t_scat /= 1000.

        # Decorrelation bandwidth (eq. 4.39)
        decorr_bw = 1.16/(2*np.pi*t_scat)
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

        t_diss, decorr_bw = go.ne2001_scint_time_bw(dist_co,
                                                    gl,
                                                    gb,
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
        snr = np.random.normal(snr, m*snr)

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
        self.s_peak_limit /= self.gain*np.sqrt(self.n_pol*self.bw*1e3)

        # Line of constant fluence
        self.fluence_limit = self.s_peak_limit / np.sqrt(w_eff)
        self.fluence_limit *= w_eff

        return self.fluence_limit
