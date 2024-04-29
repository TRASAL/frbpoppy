"""Class holding survey properties."""

import numpy as np
import os
import pandas as pd

import frbpoppy.galacticops as go
import frbpoppy.pointings as pointings
from frbpoppy.paths import paths
import frbpoppy.beam_dists as bd


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
        if self.name == 'perfect':
            self.beam_pattern = 'perfect'

        # Un-used
        self.strategy = 'regular'  # 'follow-up' not implemented
        if self.mount_type == 'transit':
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
        self.beam_size_at_fwhm = survey['beam size (deg^2)']
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

    def set_beam(self, model='perfect', size=None, random_loc=True,
                 n_sidelobes=0.5):
        """Set intensity profile.

        Set properties for int pro

        Args:
            model (str): Beam pattern. Choice from 'wsrt-apertif',
                'parkes-htru', 'chime-frb', 'gaussian', 'airy'.
            size (float): Beam size at FWHM [sq. deg].
                Defaults to that of the survey file
            random_loc (bool): Whether to calculate the precise or random
                location of each burst in the beam.
        if model in ('gaussian', 'airy'):
            n_sidelobes (int): Number of sidelobes to include. Defaults to
                cutting at the FWHM when set to 0.5.

        """
        # Set up beam properties
        self.beam_pattern = model
        self.n_sidelobes = n_sidelobes
        self.beam_array = None
        self.pixel_scale = None

        # Calculate beam properties
        if size is not None:
            self.beam_size_at_fwhm = size
        self.fwhm = 2*go.calc_sky_radius(self.beam_size_at_fwhm)

        # What should the maximum radius of the beam be?
        self.max_offset = bd.calc_max_offset(self.n_sidelobes, self.fwhm)
        self.beam_size = go.calc_sky_area(self.max_offset)

        beam_props = bd.get_beam_props(self.beam_pattern, self.fwhm)
        self.beam_size_array, self.pixel_scale, self.beam_array = beam_props

        if self.beam_size_array is not None:
            self.beam_size = self.beam_size_array
            # Not completely kosher -> a circle is not the same as a square,
            # but this is the smallest radius in which the full beam pattern
            # fits. The area of the circle extending beyond the beam pattern
            # will have an intensity of zero.
            self.max_offset = go.calc_sky_radius(self.beam_size)

        self.beam_func_oneoffs = lambda x: bd.int_pro_random(
                                 shape=x,
                                 fwhm=self.fwhm,
                                 pattern=self.beam_pattern,
                                 max_offset=self.max_offset,
                                 central_freq=self.central_freq,
                                 beam_array=self.beam_array,
                                 pixel_scale=self.pixel_scale)

        def int_pro(ra, dec, ra_p, dec_p, lst):
            return bd.int_pro_fixed(ra, dec, ra_p, dec_p, lst,
                                    pattern=self.beam_pattern,
                                    latitude=self.latitude,
                                    beam_array=self.beam_array,
                                    pixel_scale=self.pixel_scale,
                                    mount_type=self.mount_type)

        self.beam_func_rep = int_pro

    def calc_beam(self, repeaters=False, shape=None, ra=None, dec=None,
                  ra_p=None, dec_p=None, lst=None):
        """Calculate intensity profile."""
        if not repeaters and self.beam_pattern in ('airy', 'gaussian'):
            # What should the maximum radius of the beam be?
            self.max_offset = bd.calc_max_offset(self.n_sidelobes, self.fwhm)
            self.beam_size = go.calc_sky_area(self.max_offset)
        elif self.beam_pattern.startswith('perfect'):
            # What should the maximum radius of the beam be?
            self.max_offset = bd.calc_max_offset(self.n_sidelobes, self.fwhm)
            self.beam_size = go.calc_sky_area(self.max_offset)
        else:
            self.beam_size = self.beam_size_array

        if not repeaters:
            return self.beam_func_oneoffs(shape)
        else:
            return self.beam_func_rep(ra, dec, ra_p, dec_p, lst)

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
                    f_low=10e6, f_high=10e9):
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

    def calc_snr(self, s_peak, w_eff, T_sys):
        """
        Caculate the SNR of several frbs.

        Radiometer equation for single pulse (Dewey et al., 1984), but
        adapted to allow for a degradation factor reducing the peak flux
        as a pulse is stretched

        Args:
            s_peak (array): Peak flux [Jy]
            w_eff (array): Pulse width at telescope [ms]
            T_sys (array): System temperature [K]

        Returns:
            array: Signal to noise ratio based on the radiometer
                equation for a single pulse.

        """
        # Dimension check
        if s_peak.ndim == w_eff.ndim:
            pass
        elif s_peak.ndim == 1:
            s_peak = s_peak[:, None]
        elif w_eff.ndim == 1:
            w_eff = w_eff[:, None]

        if (s_peak.ndim or w_eff.ndim) > 1:
            if isinstance(T_sys, np.ndarray):
                T_sys = T_sys[:, None]

        snr = s_peak*self.gain*np.sqrt(self.n_pol*self.bw*w_eff*1e3)
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
