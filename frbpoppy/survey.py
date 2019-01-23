"""Class holding survey properties."""
import math
import numpy as np
import os
import pandas as pd
from scipy.special import j1

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
        sidelobes (int): Number of sidelobes to include
        equal_area (int/bool): Ensures beam area on sky can be equal to a beam
            pattern with a max number of sizelobes. If unwanted, set to False

    """

    def __init__(self,
                 name,
                 gain_pattern='gaussian',
                 n_sidelobes=0.5):
        """Initializing."""
        # Set up parameters
        self.name = name
        self.gain_pattern = gain_pattern
        self.n_sidelobes = n_sidelobes

        # Parse survey file
        self.read_survey_parameters()

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
        self.T_sys = survey['system temperature (K)']
        self.central_freq = int(survey['centre frequency (MHz)'])
        self.bw = survey['bandwidth (MHz)']
        self.bw_chan = survey['channel bandwidth (MHz)']
        self.n_pol = survey['number of polarizations']
        self.beam_size_fwhm = survey['beam size (deg^2)']
        self.snr_limit = survey['signal-to-noise ratio [0-1]']
        self.max_w_eff = survey['maximum pulse width (ms)']
        self.ra_min = survey['minimum RA (deg)']
        self.ra_max = survey['maximum RA (deg)']
        self.dec_min = survey['minimum DEC (deg)']
        self.dec_max = survey['maximum DEC (deg)']
        self.gl_min = survey['minimum Galactic longitude (deg)']
        self.gl_max = survey['maximum Galactic longitude (deg)']
        self.gb_min = survey['minimum Galactic latitude (deg)']
        self.gb_max = survey['maximum Galactic latitude (deg)']
        self.up_time = survey['fractional uptime [0-1]']

    def in_region(self, frbs):
        """
        Check if the given frbs are within the survey region.

        Args:
            frbs (Frbs): Frbs of which to check whether in survey region

        Returns:
            array: Boolean mask denoting whether frbs are within survey region

        """
        # Create mask with False
        mask = np.ones_like(frbs.ra, dtype=bool)

        # Ensure in correct format
        frbs.gl[frbs.gl > 180.] -= 360.

        # Create region masks
        gl_limits = (frbs.gl > self.gl_max) | (frbs.gl < self.gl_min)
        gb_limits = (frbs.gb > self.gb_max) | (frbs.gb < self.gb_min)
        ra_limits = (frbs.ra > self.ra_max) | (frbs.ra < self.ra_min)
        dec_limits = (frbs.dec > self.dec_max) | (frbs.dec < self.dec_min)
        mask[gl_limits] = False
        mask[gb_limits] = False
        mask[ra_limits] = False
        mask[dec_limits] = False

        return mask

    def max_offset(self, x):
        """Calculate the maximum offset of an FRB in an Airy disk.

        Args:
            x (int): Maximum sidelobe wanted

        """
        # Null points of kasin for allow a number of sidelobes
        kasin_nulls = [3.831706,
                       7.015587,
                       10.173468,
                       13.323692,
                       16.47063,
                       19.615859,
                       22.760084,
                       25.903672,
                       29.046829,
                       32.18968,
                       35.332308,
                       38.474766]

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

    def intensity_profile(self, n_gen=1, dimensions=2):
        """Calculate intensity profile."""
        # Calculate Full Width Half Maximum from beamsize
        self.fwhm = 2*math.sqrt(self.beam_size_fwhm/math.pi) * 60  # [arcmin]
        offset = self.fwhm/2  # Radius = diameter/2.

        if dimensions == 2:  # 2D
            offset *= np.sqrt(np.random.random(n_gen))
        elif dimensions == 1:  # 1D
            offset *= np.random.random(n_gen)

        # Allow for a perfect beam pattern in which all is detected
        if self.gain_pattern == 'perfect':
            int_pro = np.ones(n_gen)
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
            shape = beam_array.shape
            ran_x = np.random.randint(0, shape[0], n_gen)
            ran_y = np.random.randint(0, shape[1], n_gen)
            int_pro = beam_array[ran_x, ran_y]
            offset = np.sqrt((ran_x-shape[0]/2)**2 + (ran_y-shape[1]/2)**2)

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

    def dm_smear(self, frbs, dm_err=0.2):
        """
        Calculate delay in pulse across a channel due to dm smearing.

        Formula's based on 'Handbook of Pulsar Astronomy" by Duncan Lorimer
        & Michael Kramer, section A2.4. Note the power of the forefactor has
        changed due to the central frequency being given in MHz.

        Args:
            frbs (FRBs): FRB object with a dm attribute
            dm_err (float): Error on dispersion measure. Defaults to 0.2

        Returns:
            t_dm, t_dm_err (arrays): Time of delay [ms] at central band
                frequency, with its error assuming a 20% uncertainty in the
                dispersion measure

        """
        t_dm = 8.297616e6 * self.bw_chan * frbs.dm * (self.central_freq)**-3
        t_dm_err = t_dm*dm_err

        return t_dm, t_dm_err

    def calc_scat(self, dm):
        """Calculate scattering timescale for FRBs.

        Offset according to Lorimer et al. (doi:10.1093/mnrasl/slt098)

        Args:
            dm (array): Dispersion Measure

        Returns:
            array: Scattering timescales [ms]

        """
        freq = self.central_freq
        t_scat = go.scatter_bhat(dm, scindex=-3.86, offset=-9.5, freq=freq)
        return t_scat

    def calc_Ts(self, frbs):
        """Set temperatures for frbs."""
        T_sky = self.calc_T_sky(frbs)
        T_tot = self.T_sys + T_sky
        return T_sky, T_tot

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
        T_sky_haslam = np.take(T_sky_list, index)

        # scale temperature
        # Assuming dominated by syncrotron radiation
        T_sky = T_sky_haslam * (self.central_freq/408.0)**(-2.6)

        return T_sky

    def calc_s_peak(self, frbs, f_low=10e6, f_high=10e9):
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
        sp = frbs.si + 1
        sm = frbs.si - 1

        # Convert distance in Gpc to metres
        dist = frbs.dist_co * 1e9 * 3.0856775814913673 * 1e16

        # Convert luminosity to Watts
        lum = frbs.lum_bol * 1e-7

        freq_frac = (f_2**sp - f_1**sp) / (f_2 - f_1)
        nom = lum * (1+frbs.z)**sm * freq_frac
        den = 4*np.pi*dist**2 * (f_high**sp - f_low**sp)
        s_peak = nom/den

        # Convert to Janskys
        s_peak *= 1e26

        return s_peak

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
        w_eff = np.sqrt(frbs.w_arr**2 +
                        frbs.t_dm**2 +
                        frbs.t_dm_err**2 +
                        frbs.t_scat**2 +
                        self.t_samp**2)
        return w_eff

    def calc_snr(self, frbs):
        """
        Caculate the SNR of several frbs.

        Args:
            frbs (FRBs): FRBs of which to calculate the signal to noise

        Returns:
            array: Signal to noise ratio based on the radiometer
                equation for a single pulse.

        """
        # Radiometer equation for single pulse (Dewey et al., 1984)
        sp = frbs.s_peak
        snr = sp*self.gain*np.sqrt(self.n_pol*self.bw*frbs.w_eff*1e3)
        snr /= (frbs.T_tot * self.beta)
        return snr

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
        self.s_peak_limit = self.snr_limit*self.T_sys*self.beta
        self.s_peak_limit /= self.gain*math.sqrt(self.n_pol*self.bw*1e3)

        # Line of constant fluence
        self.fluence_limit = self.s_peak_limit / math.sqrt(w_eff)
        self.fluence_limit *= w_eff

        return self.fluence_limit
