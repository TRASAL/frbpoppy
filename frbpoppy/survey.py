import math
import os
import random
from scipy.special import j1

import frbpoppy.galacticops as go
from frbpoppy.log import pprint
from frbpoppy.paths import paths


class Rates:
    """
    Class to hold rate counters

    Args:
        det (int, optional): Number detected
        faint (int, optional): Number to faint to detect
        jy (int, optional): Number arrived within survey > 1 Jy
        out (int, optional): Number outside survey, space or timewise
        sky (int, optional): Number per sky above > 1 Jy
        vol (int, optional): Number per Gpc^3

    """

    def __init__(self):

        # Rates
        self.det = 0
        self.faint = 0
        self.jy = 0
        self.out = 0
        self.sky = 0
        self.vol = 0

    def tot(self):
        """Calculate the total number of rates"""
        return self.det + self.out + self.faint


class Survey:
    """
    Method containing survey parameters and functions

    Args:
        survey_name (str): Name of survey with which to observe population. Can
            either be a predefined survey present in frbpoppy
            (see data/surveys/) or a path name to a new survey
            filename
        pattern (str): Set gain pattern to be either 'gaussian' or 'airy'.
            Defaults to 'gaussian'
    """

    def __init__(self, survey_name, pattern='gaussian'):

        # Find survey file
        if os.path.isfile(survey_name):
            f = open(survey_name, 'r')
        else:
            # Find standard survey files
            try:
                path = os.path.join(paths.surveys(), survey_name)
                f = open(path, 'r')
            except IOError:
                s = 'Survey file {0} does not exist'.format(survey_name)
                raise IOError(s)

        # Parse survey file
        self.gl2_min = None
        self.gl2_max = None
        self.parse(f)
        self.discoveries = 0
        self.survey_name = survey_name
        self.pointings_list = None
        self.gains_list = None
        self.t_obs_list = None
        self.T_sky_list = go.load_T_sky()
        self.gain_pattern = pattern
        self.aa = False  # Whether aperture array

        # Counters
        self.frb_rates = Rates()
        self.src_rates = Rates()

    def __str__(self):
        """Define how to print a survey object to a console"""

        s = 'Survey properties:'

        attributes = []
        for e in self.__dict__:
            attr = '\n\t{0:13.12}{1:.60}'.format(e, str(self.__dict__[e]))
            attributes.append(attr)

        s += ''.join(attributes)

        return s

    def parse(self, f):
        """
        Attempt to parse an already opened survey file

        Args:
            f (str): Filename, see Survey class

        Returns:
            Various attributes
        """

        for line in f:

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
                self.beam_size = float(v)  # [deg**2]
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
            elif p.count('uptime'):
                self.uptime = float(v)  # % of time in use
                if self.uptime > 1.0:
                    self.uptime = 1.0
            elif p.count('signal-to-noise'):
                self.snr_limit = float(v)  # Minimum snr required for detection
            elif p.count('gain pattern'):
                self.gain_pat = v  # Gain pattern of telescope
            elif p.count('Aperture Array'):
                self.aa = True
            elif p.count('reference'):
                pass
            elif p.count('colour'):
                pass
            else:
                pprint('Parameter {0} not recognised'.format(p))

        f.close()

    def in_region(self, src):
        """
        Check if a given source is within the survey region

        Args:
            src (Source): Source of which to check whether in survey region

        Returns:
            True: If source is within survey region
            False: If source is outside survey region
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

    def intensity_profile(self):
        """Calculate intensity profile"""

        # Angular variable on sky, defined so that at fwhm/2, the
        # intensity profile is exactly 0.5. It's multiplied by the sqrt of
        # a random number to ensure the distribution of variables on the
        # sky remains uniform, not that it increases towards the centre
        # (uniform random point within circle). You could see this as
        # the offset from the centre of the beam.

        self.fwhm = 2*math.sqrt(self.beam_size/math.pi) * 60  # [arcmin]

        offset = self.fwhm * math.sqrt(random.random()) / 2.0

        if self.gain_pattern == 'gaussian':

            # Formula's based on 'Interferometry and Synthesis in Radio
            # Astronomy' by A. Richard Thompson, James. M. Moran and
            # George W. Swenson, JR. (Second edition), around p. 15

            # Intensity profile
            alpha = 2*math.sqrt(math.log(2))
            int_pro = math.exp(-(alpha*offset/self.fwhm)**2)

        if self.gain_pattern == 'airy':
            c = 299792458
            conv = math.pi/(60*180.)  # Conversion arcmins -> radians
            eff_diam = c/(self.central_freq*1e6*conv*self.fwhm)
            a = eff_diam/2  # Effective radius of telescope
            lamda = c/(self.central_freq*1e6)
            kasin = (2*math.pi*a/lamda)*math.sin(offset*conv)
            int_pro = 4*(j1(kasin)/kasin)**2

        if self.gain_pattern == 'tophat':
            int_pro = 1
            if random.random() > 0.5:
                int_pro = 0

        return int_pro

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
                frequency, with its error assuming a
                20% uncertainty in the dispersion measure

        """
        t_dm = 8.297616e6 * self.bw_chan * src.dm * (self.central_freq)**-3
        t_dm_err = (t_dm/src.dm)*(dm_err*src.dm)

        src.t_dm = t_dm
        src.t_dm_err = t_dm_err

    def scat(self, src):
        """
        Set scattering timescale for source
        """
        # Offset according to Lorimer et al. (doi:10.1093/mnrasl/slt098)
        src.t_scat = go.scatter_bhat(src.dm,
                                     scindex=-3.86,
                                     offset=-9.5,
                                     freq=self.central_freq)

    def calc_Ts(self, src):
        """Set temperatures for source"""
        self.calc_T_sky(src)
        src.T_tot = self.T_sys + src.T_sky

    def calc_T_sky(self, src):
        """
        Calculate the sky temperature from the Haslam table, before scaling to
        the survey frequency. The temperature sky map is given in the weird
        units of HealPix and despite looking up info on this coordinate system,
        I don't have the foggiest idea of how to transform these to galactic
        coordinates. I have therefore directly copied the following code from
        psrpoppy in the assumption Sam Bates managed to figure it out.

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
            f_low (float): Source emission lower frequency limit [Hz]. Defaults
                to 10e6
            f_high (float): Source emission higher frequency limit [Hz].
                Defaults to 10e6

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

        # Convert distance to metres
        dist = src.dist * 3.08567758149137e25

        # Convert luminosity to Watts
        lum = frb.lum_bol * 1e-7

        freq_frac = (f_2**sp - f_1**sp) / (f_2 - f_1)
        nom = lum * (1+src.z)**sm * freq_frac
        den = 4*math.pi*dist**2 * (f_high**sp - f_low**sp)
        s_peak = nom/den

        # Convert to Janskys
        s_peak *= 1e26

        frb.s_peak = s_peak

    def obs_prop(self, frb, src, pop):
        """
        Set various observation properties of an FRB

        Args:
            frb (class): FRB of which to calculate the signal to noise
            src (class): Source to which the FRB belongs
            pop (class): Population to which the source belongs

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
        frb.w_eff = math.sqrt(frb.w_int**2 + src.t_dm**2 + src.t_dm_err**2 +
                              src.t_scat**2 + self.t_samp**2)

        # Calculate flux density
        self.calc_s_peak(frb, src, f_low=pop.f_min, f_high=pop.f_max)

        # Radiometer equation for single pulse (Dewey et al., 1984)
        sp = frb.s_peak
        frb.snr = sp*self.gain*math.sqrt(self.n_pol*self.bw*frb.w_eff*1e3)
        frb.snr /= (src.T_tot * self.beta)

        # Calculate fluence [Jy*ms]
        frb.fluence = frb.s_peak * frb.w_eff

        # Account for offset in beam
        frb.snr *= self.intensity_profile()

    def scint(self, frb, src):
        """
        Calculate scintillation effect on the signal to noise ratio (rather
        than adapting the flux, as the snr can change per survey attempt).
        Formulas based on 'Handbook of Pulsar Astronomy" by Duncan Lorimer &
        Michael Kramer, section 4.2.

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

            t_diss, decorr_bw = go.ne2001_scint_time_bw(src.dist,
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

    def scale_rates(self, pop):
        """
        Scale all rates according to integration time and uptime.

        Args:
            pop (Population): Population class
        """
        rates = (self.frb_rates, self.src_rates)
        n = 0

        for r in rates:
            area_sky = 4*math.pi*(180/math.pi)**2   # In sq. degrees
            f_area = (self.beam_size * r.tot()) / ((r.det + r.faint)*area_sky)
            f_time = 86400 / self.t_obs  # pop.time
            det = r.det * f_area * f_time
            faint = r.faint * f_area * f_time
            out = r.out + r.det - det + r.faint - faint

            vol = r.tot() / pop.v_max * (365.25*86400/pop.time)

            s_rates = Rates()
            s_rates.det = det
            s_rates.faint = faint
            s_rates.out = out
            s_rates.vol = vol

            if n == 0:
                self.s_frb_rates = s_rates
                n += 1
            else:
                self.s_src_rates = s_rates

    def rates(self, pop, output=True, scaled=True):
        """
        Print survey rates, such as the number of detected sources etc.

        Args:
            output (bool, optional): Whether to print out the rates or not
            scaled (bool, optional): Print scaled (default) or normal rates
        """
        f = self.s_frb_rates
        s = self.s_src_rates
        if not scaled:
            f = self.frb_rates
            s = self.src_rates

        r = '{:20.19} {:>10} {:>10} {:>10}\n'

        # Set up title
        days = (pop.time/86400)

        t = r.format(self.survey_name, 'Days', 'FRBs', 'Sources')
        line = '-'*len(t.split('\n')[-2].strip()) + '\n'
        t += line

        tot = ('In population', round(days), round(f.tot()), round(s.tot()))
        det = ('Detected', round(days), round(f.det), round(s.det))
        faint = ('Too faint', round(days), round(f.faint), round(s.faint))
        out = ('Outside survey', round(days), round(f.out), round(s.out))
        vol = ('/Gpc^3', 365.25, round(f.vol), round(s.vol))
        if f.det > 0:
            exp = round(days/f.det, 4)
        else:
            exp = '?'
        exp_frb = ('Expected', exp, 1, '-')

        t += r.format(*tot)
        t += r.format(*det)
        t += r.format(*faint)
        t += r.format(*out)
        t += r.format(*vol)
        t += r.format(*exp_frb)
        t += line

        for n in t.split('\n')[:-1]:
            if output:
                pprint(n)
