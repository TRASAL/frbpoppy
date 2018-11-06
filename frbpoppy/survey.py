"""Class holding survey properties."""
import math
import numpy as np
import os
import random
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
                 sidelobes=1,
                 equal_area=False):
        """Initializing."""
        # Set up parameters
        self.name = name
        self.gain_pattern = gain_pattern
        self.sidelobes = sidelobes
        self.equal_area = equal_area

        self.beam_array = None
        self.T_sky_list = go.load_T_sky()

        # Parse survey file
        self.find_survey_file()
        self.parse()

        # Beam pattern arguments
        self._kasin_nulls = [3.831706,
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

        # Set beam file so that it is only imported once
        if gain_pattern in ['parkes', 'apertif']:
            place = paths.models() + f'/beams/{gain_pattern}.npy'
            self.beam_array = np.load(place)

    def __str__(self):
        """Define how to print a survey object to a console."""
        s = 'Survey properties:'

        attributes = []
        for e in self.__dict__:
            attr = '\n\t{0:13.12}{1:.60}'.format(e, str(self.__dict__[e]))
            attributes.append(attr)

        s += ''.join(attributes)

        return s

    def find_survey_file(self):
        """Use defaults to find survey file."""
        if os.path.isfile(self.name):
            path = self.name
        else:
            # Find standard survey files
            try:
                path = os.path.join(paths.surveys(), self.name)
            except IOError:
                s = 'Survey file {0} does not exist'.format(self.name)
                raise IOError(s)

        self.survey_file = path

    def parse(self):
        """Attempt to parse an already opened survey file."""
        open_file = open(self.survey_file, 'r')
        for line in open_file:

            # Ignore comments
            if line[0] == '#':
                continue

            # Parse arguments
            a = line.split('!')
            p = a[1].strip()
            v = a[0].strip()

            if p.count('survey degradation'):
                self.beta = float(v)
            elif p.count('antenna gain'):
                self.gain = float(v)  # [K/Jy]
            elif p.count('integration time'):
                self.t_obs = float(v)  # [s]
            elif p.count('sampling'):
                self.t_samp = float(v)  # [ms]
            elif p.count('system temperature'):
                self.T_sys = float(v)  # [K]
            elif p.count('centre frequency'):
                self.central_freq = float(v)  # [MHz]
            elif p.startswith('bandwidth'):
                self.bw = float(v)  # [MHz]
            elif p.count('channel bandwidth'):
                self.bw_chan = float(v)  # [MHz]
            elif p.count('polarizations'):
                self.n_pol = float(v)  # number of polarizations
            elif p.count('beam size'):
                self.beam_size_fwhm = float(v)  # [deg**2]
            elif p.count('maximum pulse width'):
                self.max_w_eff = float(v)  # [ms]
            elif p.count('minimum RA'):
                self.ra_min = float(v)  # [deg]
            elif p.count('maximum RA'):
                self.ra_max = float(v)  # [deg]
            elif p.count('minimum DEC'):
                self.dec_min = float(v)  # [deg]
            elif p.count('maximum DEC'):
                self.dec_max = float(v)  # [deg]
            elif p.count('minimum Galactic longitude'):
                self.gl_min = float(v)  # [deg]
            elif p.count('maximum Galactic longitude'):
                self.gl_max = float(v)  # [deg]
            elif p.count('minimum 2 Galactic longitude'):
                self.gl2_min = float(v)  # [deg]
            elif p.count('maximum 2 Galactic longitude'):
                self.gl2_max = float(v)  # [deg]
            elif p.count('minimum Galactic latitude'):
                self.gb_min = float(v)  # [deg]
            elif p.count('maximum Galactic latitude'):
                self.gb_max = float(v)  # [deg]
            elif p.count('signal-to-noise'):
                self.snr_limit = float(v)  # Minimum s/r required for detection
            elif p.count('gain pattern'):
                pass  # Gain pattern of telescope
            elif p.count('Aperture Array'):
                pass
            elif p.count('uptime'):
                pass
            elif p.count('reference'):
                pass
            elif p.count('colour'):
                pass
            else:
                pprint('Parameter {0} not recognised'.format(p))

        open_file.close()

    def in_region(self, src):
        """
        Check if a given source is within the survey region.

        Args:
            src (Source): Source of which to check whether in survey region

        Returns:
            bool: If source is within a survey region

        """
        if src.gl > 180.:
            src.gl -= 360.

        if src.gl > self.gl_max or src.gl < self.gl_min:
            if self.gl2_min:
                if src.gl > self.gl2_max or src.gl < self.gl2_min:
                    return False
            else:
                return False
        if src.gb > self.gb_max or src.gb < self.gb_min:
            return False

        if src.ra > self.ra_max or src.ra < self.ra_min:
            return False
        if src.dec > self.dec_max or src.dec < self.dec_min:
            return False

        return True

    def max_offset(self, x):
        """Calculate the maximum offset of an FRB in an Airy disk.

        Args:
            x (int): Maximum sidelobe wanted

        """
        # Allow for cut at FWHM
        if x == 0.5:
            return 1

        try:
            arcsin = math.asin(self.fwhm*self._kasin_nulls[x]/(60*180))
        except ValueError:
            m = f'Beamsize including sidelobes would be larger than sky \n'
            A = (90/self._kasin_nulls[x])**2*math.pi
            m += f'Ensure beamsize is smaller than {A}'
            raise ValueError(m)

        return 2/self.fwhm * 60*180/math.pi * arcsin

    def intensity_profile(self, sidelobes=1, test=False, equal_area=False,
                          dimensions=2):
        """Calculate intensity profile."""
        self.fwhm = 2*math.sqrt(self.beam_size_fwhm/math.pi) * 60  # [arcmin]

        offset = self.fwhm/2  # Radius = Diameter/2.

        if self.gain_pattern == 'perfect':
            int_pro = 1
            self.beam_size = self.beam_size_fwhm
        else:
            max_offset = offset*self.max_offset(sidelobes)
            self.beam_size = math.pi*(max_offset/60)**2  # [sq degrees]

        if dimensions == 2:  # 2D
            offset *= math.sqrt(random.random())
        elif dimensions == 1:  # 1D
            offset *= random.random()

        if self.gain_pattern == 'tophat':
            int_pro = 1
            if random.random() > 0.5:
                int_pro = 0

        # Formula's based on 'Interferometry and Synthesis in Radio
        # Astronomy' by A. Richard Thompson, James. M. Moran and
        # George W. Swenson, JR. (Second edition), around p. 15

        elif self.gain_pattern == 'gaussian':
            # Set the maximum offset equal to the null after a sidelobe
            # I realise this pattern isn't an airy, but you have to cut
            # somewhere
            offset *= self.max_offset(sidelobes)

            alpha = 2*math.sqrt(math.log(2))
            int_pro = math.exp(-(alpha*offset/self.fwhm)**2)

        elif self.gain_pattern == 'airy':
            int_pro = None

            if type(equal_area) == int:
                offset *= self.max_offset(equal_area)
                # Check whether offset is within the sidelobe you're including
                if offset > self.max_offset(sidelobes)*self.fwhm/2.:
                    int_pro = 0
                # Calculate beamsize
                max_offset = self.max_offset(equal_area)*self.fwhm/2.
                self.beam_size = math.pi*(max_offset/60)**2  # [sq degrees]
            else:
                # Set the maximum offset equal to the null after a sidelobe
                offset *= self.max_offset(sidelobes)

            if int_pro is None:
                c = 299792458
                conv = math.pi/(60*180)  # Conversion arcmins -> radians
                eff_diam = c/(self.central_freq*1e6*conv*self.fwhm)
                a = eff_diam/2  # Effective radius of telescope
                lamda = c/(self.central_freq*1e6)
                ka = (2*math.pi*a/lamda)
                kasin = ka*math.sin(offset*conv)
                int_pro = 4*(j1(kasin)/kasin)**2

        elif self.gain_pattern in ['parkes', 'apertif']:
            shape = self.beam_array.shape
            ran_x = np.random.randint(0, shape[0])
            ran_y = np.random.randint(0, shape[1])
            int_pro = self.beam_array[ran_x, ran_y]
            offset = math.sqrt((ran_x-shape[0]/2)**2 + (ran_y-shape[1]/2)**2)

            # Scaling factors to correct for pixel scale
            if self.gain_pattern == 'apertif':  # 1 pixel = 0.94'
                offset *= 240/256  # [arcmin]
                self.beam_size = 25.
            if self.gain_pattern == 'parkes':  # 1 pixel = 54"
                offset *= 0.9  # [arcmin]
                self.beam_size = 9.

        if test is False:
            return int_pro
        else:
            return offset, int_pro

    def dm_smear(self, src, dm_err=0.2):
        """
        Calculate delay in pulse across a channel due to dm smearing.

        Formula's based on 'Handbook of Pulsar Astronomy" by Duncan Lorimer
        & Michael Kramer, section A2.4. Note the power of the forefactor has
        changed due to the central frequency being given in MHz.

        Args:
            src (Source): Source object with a dm attribute
            dm_err (float): Error on dispersion measure. Defaults to 0.2

        Returns:
            t_dm, t_dm_err (float): Time of delay [ms] at central band
                frequency, with its error assuming a 20% uncertainty in the
                dispersion measure

        """
        t_dm = 8.297616e6 * self.bw_chan * src.dm * (self.central_freq)**-3
        t_dm_err = t_dm*dm_err

        src.t_dm = t_dm
        src.t_dm_err = t_dm_err

    def scat(self, src):
        """Set scattering timescale for source."""
        # Offset according to Lorimer et al. (doi:10.1093/mnrasl/slt098)
        src.t_scat = go.scatter_bhat(src.dm,
                                     scindex=-3.86,
                                     offset=-9.5,
                                     freq=self.central_freq)

    def calc_Ts(self, src):
        """Set temperatures for source."""
        self.calc_T_sky(src)
        src.T_tot = self.T_sys + src.T_sky

    def calc_T_sky(self, src):
        """
        Calculate the sky temperature from the Haslam table.

        Afterwards scale to the survey frequency. The temperature sky map is
        given in the weird units of HealPix and despite looking up info on this
        coordinate system, I don't have the foggiest idea of how to transform
        these to galactic coordinates. I have therefore directly copied the
        following code from psrpoppy in the assumption Sam Bates managed to
        figure it out.

        Args:
            src (class): Needed for coordinates
        Returns:
            src.T_sky (float): Sky temperature [K]
        """
        # ensure l is in range 0 -> 360
        B = src.gb
        if src.gl < 0.:
            L = 360 + src.gl
        else:
            L = src.gl

        # convert from l and b to list indices
        j = B + 90.5
        if j > 179:
            j = 179

        nl = L - 0.5
        if L < 0.5:
            nl = 359
        i = float(nl) / 4.

        T_sky_haslam = self.T_sky_list[180*int(i) + int(j)]

        # scale temperature
        # Assuming dominated by syncrotron radiation
        T_sky = T_sky_haslam * (self.central_freq/408.0)**(-2.6)

        src.T_sky = T_sky

    def calc_s_peak(self, frb, src, f_low=10e6, f_high=10e9):
        """
        Calculate the mean spectral flux density.

        Following Lorimer et al, 2013, eq. 9., at the central frequency
        of the survey.

        Args:
            frb (class): FRB
            src (class): Source of FRB
            f_low (float): Source emission lower frequency limit [Hz].
            f_high (float): Source emission higher frequency limit [Hz].

        Returns:
            frb.s_peak (float): Mean spectral flux density [Jy]

        """
        # Limits observing bandwidth (as seen in rest frame source)
        f_1 = (self.central_freq - 0.5*self.bw)
        f_1 *= 1e6  # MHz -> Hz
        f_2 = (self.central_freq + 0.5*self.bw)
        f_2 *= 1e6  # MHz -> Hz

        # Spectral index
        sp = frb.si + 1
        sm = frb.si - 1

        # Convert distance in Gpc to metres
        dist = src.dist_co * 1e9 * 3.0856775814913673 * 1e16

        # Convert luminosity to Watts
        lum = frb.lum_bol * 1e-7

        freq_frac = (f_2**sp - f_1**sp) / (f_2 - f_1)
        nom = lum * (1+src.z)**sm * freq_frac
        den = 4*math.pi*dist**2 * (f_high**sp - f_low**sp)
        s_peak = nom/den

        # Convert to Janskys
        s_peak *= 1e26

        frb.s_peak = s_peak

    def obs_prop(self, frb, src, f_min, f_max):
        """
        Set various observation properties of an FRB.

        Args:
            frb (class): FRB of which to calculate the signal to noise
            src (class): Source to which the FRB belongs
            f_min (float): Source emission lower frequency limit [Hz].
            f_max (float): Source emission higher frequency limit [Hz].

        Returns:
            frb.snr (float): Signal to noise ratio based on the radiometer
                equation for a single pulse. Will return -2.0 if src is not
                in survey region
            frb.w_eff (float): Observed pulse width [ms]. Will return 0. if
                source not in survey region
            frb.s_peak (float): Mean spectral flux density per observed
                source [Jy]
            frb.fluence (float): Fluence of the observed pulse [Jy*ms]

        """
        # Effective pulse width [ms]
        # From Narayan (1987, DOI: 10.1086/165442)
        # Also Cordes & McLaughlin (2003, DOI: 10.1086/378231)
        # For details see p. 30 of Emily Petroff's thesis (2016), found here:
        # http://hdl.handle.net/1959.3/417307
        if not frb.w_eff:  # Hack to go from frbcat df to pop
            frb.w_eff = math.sqrt(frb.w_arr**2 + src.t_dm**2 +
                                  src.t_dm_err**2 + src.t_scat**2 +
                                  self.t_samp**2)

        # Calculate flux density
        if not frb.s_peak:  # Hack to go from frbcat df to pop
            self.calc_s_peak(frb, src, f_low=f_min, f_high=f_max)

        # Account for offset in beam
        frb.s_peak *= self.intensity_profile(sidelobes=self.sidelobes,
                                             equal_area=self.equal_area)

        # Radiometer equation for single pulse (Dewey et al., 1984)
        sp = frb.s_peak
        frb.snr = sp*self.gain*math.sqrt(self.n_pol*self.bw*frb.w_eff*1e3)
        frb.snr /= (src.T_tot * self.beta)

        # Calculate fluence [Jy*ms]
        frb.fluence = frb.s_peak * frb.w_eff

    def scint(self, frb, src):
        """
        Calculate scintillation effect on the signal to noise ratio.

        (Rather than adapting the flux, as the snr can change per survey
        attempt). Formulas based on 'Handbook of Pulsar Astronomy" by Duncan
        Lorimer & Michael Kramer, section 4.2.

        Args:
            frb (class): FRB
            src (class): Source of FRB
        Returns:
            frb.snr (float): Signal to noise ratio modulated by scintillation

        """
        # Calculate scattering
        self.scat(src)
        # Convert to seconds
        t_scat = src.t_scat
        t_scat /= 1000.

        # Decorrelation bandwidth (eq. 4.39)
        decorr_bw = 1.16/(2*math.pi*t_scat)
        # Convert to MHz
        decorr_bw /= 1e6

        # Scintillation strength (eq. 4.33)
        u = math.sqrt(self.central_freq / decorr_bw)

        # Strong scintillation
        if u < 1:
            # (eq. 4.35)
            m = math.sqrt(u**(5/3))

        # Weak scintillation
        else:

            # Refractive scintillation (eq. 4.47)
            m_riss = u**-(1/3)

            # Taking the average kappa value
            kappa = 0.15

            t_diss, decorr_bw = go.ne2001_scint_time_bw(src.dist_co,
                                                        src.gl,
                                                        src.gb,
                                                        self.central_freq)

            # Following Cordes and Lazio (1991) (eq. 4.43)
            if t_diss is None:
                n_t = 1.
            else:
                n_t = 1 + kappa * self.t_obs / t_diss

            if decorr_bw is None:
                n_f = 1.
            else:
                n_f = 1 + kappa * self.bw / decorr_bw

            # Diffractive scintillation (eq. 4.41)
            m_diss = 1 / math.sqrt(n_t * n_f)

            # (eq. 4.48)
            m = math.sqrt(m_diss**2 + m_riss**2 + m_diss*m_riss)

        # Distribute the scintillation according to gaussian distribution
        frb.snr = random.gauss(frb.snr, m*frb.snr)

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
