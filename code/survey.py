import math
import os
import random

import galacticops as go
from log import pprint


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
        self.det = 0
        self.faint = 0
        self.jy = 0
        self.out = 0
        self.sky = 0
        self.vol = 0

        self.s_det = 0
        self.s_faint = 0
        self.s_jy = 0
        self.s_out = 0
        self.s_sky = 0
        self.s_vol = 0

    def tot(self):
        """Calculate the total number of rates"""
        return self.det + self.out + self.faint

    def s_tot(self):
        """Calculate the total number of scaled rates"""
        return self.s_det + self.s_out + self.s_faint

class Survey:
    """
    Method containing survey parameters and functions

    Args:
        survey_name (str): Name of survey with which to observe population. Can
            either be a predefined survey present in frbpoppy
            (see data/surveys/) or a path name to a new survey
            filename
        pattern (str): Set gain pattern. Defaults to 'gaussian'
    """

    def __init__(self, survey_name, pattern='gaussian'):

        # Find survey file
        if os.path.isfile(survey_name):
            f = open(survey_name, 'r')
        else:
            # Find standard survey files
            try:
                survey_dir = os.path.dirname(__file__) + '/../data/surveys/'
                path = os.path.join(survey_dir, survey_name)
                f = open(path, 'r')
            except IOError:
                s = 'Survey file {0} does not exist'.format(survey_name)
                raise IOError(s)

        # Parse survey file
        self.parse(f)
        self.discoveries = 0
        self.survey_name = survey_name
        self.pointings_list = None
        self.gains_list = None
        self.t_obs_list = None
        self.T_sky_list = go.load_T_sky()
        self.gain_pattern = pattern
        self.aa = False  # Whether aperture array
        self.int_pro = self.intensity_profile()

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
            elif p.count('minimum Galactic'):
                self.gl_min = float(v)  # [deg]
            elif p.count('maximum Galactic'):
                self.gl_max = float(v)  # [deg]
            elif p.count('minimum abs'):
                self.gb_min = float(v)  # min(abs(galactic latitude)) [deg]
            elif p.count('maximum abs'):
                self.gb_max = float(v)  # max(abs(galactic latitude)) [deg]
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
            return False

        abs_gb = math.fabs(src.gb)
        if abs_gb > self.gb_max or abs_gb < self.gb_min:
            return False

        ra, dec = go.lb_to_radec(src.gl, src.gb)

        if ra > self.ra_max or ra < self.ra_min:
            return False
        if dec > self.dec_max or dec < self.dec_min:
            return False

        return True

    def intensity_profile(self):
        """Calculate intensity profile"""

        if self.gain_pattern == 'gaussian':

            # Formula's based on 'Interferometry and Synthesis in Radio
            # Astronomy' by A. Richard Thompson, James. M. Moran and
            # George W. Swenson, JR. (Second edition), around p. 15

            # Angular variable on sky, defined so that at fwhm/2, the
            # intensity profile is exactly 0.5. It's multiplied by the sqrt of
            # a random number to ensure the distribution of variables on the
            # sky remains uniform, not that it increases towards the centre
            # (uniform random point within circle). You could see this as
            # the calculated offset from the centre of the beam. The fwhm
            # has been dropped as it is also present in the intensity profile
            xi = math.sqrt(random.random()) / 2.

            # Intensity profile
            alpha = 2*math.sqrt(math.log(2))
            int_pro = math.exp(-(alpha*xi)**2)

        return int_pro

    def dm_smear(self, src, dm_err=0.2):
        """
        Calculate delay in pulse across a channel due to dm smearing. Formula's
        based on 'Handbook of Pulsar Astronomy" by Duncan Lorimer & Michael
        Kramer, section A2.4. Note the power of the forefactor has changed due
        to the central frequency being given in MHz.

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
        Calculate the mean spectral flux density following Lorimer et al, 2013,
        eq. 9., at the central frequency of the survey.

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
        frb.snr *= self.int_pro

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
        Scale all rates according to intergration time and uptime

        Args:
            pop (Population): Population class
        """

        for r in (self.frb_rates, self.src_rates):

            f = (self.t_obs/pop.time)
            det = r.det * self.uptime * f
            faint = r.faint * self.uptime * f
            out = r.out + r.det - det + r.faint - faint

            vol = r.tot() / pop.v_max / (pop.time/86400)

            r.s_det = round(det)
            r.s_faint = round(faint)
            r.s_out = round(out)
            r.s_vol = round(vol)

    def rates(self, pop):
        """
        Print survey rates, such as the number of detected sources etc

        Args:
            pop (Population): Population class
        """

        t = self.survey_name + ':\n'

        f = self.frb_rates
        s = self.src_rates

        r = '   {:15.14} {:>10} {:>10}\n'
        t += r.format('Number of', 'FRBs', 'Sources')
        t += '   ' + '-'*len(t.split('\n')[-2].strip()) + '\n'
        t += r.format('in population', f.s_tot(), s.s_tot())
        t += r.format('detected', f.s_det, s.s_det)
        t += r.format('too faint', f.s_faint, s.s_faint)
        t += r.format('outside survey', f.s_out, s.s_out)
        t += r.format('/Gpc^3/day', f.s_vol, s.s_vol)

        for n in t.split('\n')[:-1]:
            pprint(n)