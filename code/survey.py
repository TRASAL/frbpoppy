import math
import os
import random

import galacticops as go


class Survey:
    """
    Method containing survey parameters and functions

    Args:
        survey_name (str): Name of survey with which to observe population. Can
                           either be a predefined survey present in frbpoppy
                           (see data/surveys/) or a path name to a new survey
                           filename
        pattern (str): Set gain pattern
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

        # Counters
        self.n_det = 0  # Number of detected sources
        self.n_faint = 0  # Number of sources too faint to detect
        self.n_out = 0  # Number of sources outside detection region

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
        Attempt to parse survey file already opened

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
            elif p.count('half maximum'):
                self.fwhm = float(v)  # [arcmin]
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
            elif p.count('coverage'):
                self.coverage = float(v)  # coverage % of sky survey area
                if self.coverage > 1.0:
                    self.coverage = 1.0
            elif p.count('signal-to-noise'):
                self.snr_limit = float(v)
            elif p.count('gain pattern'):
                self.gain_pat = v  # default = gaussian
            elif p.count('Aperture Array'):
                self.aa = True
            else:
                print('Parameter {0} not recognised'.format(p))

        f.close()

    def in_region(self, source):
        """
        Check if a given source is within the survey region

        Args:
            source (class): Source of which to check whether in survey region

        Returns:
            True: If source is within survey region
            False: If source is outside survey region
        """

        if source.gl > 180.:
            source.gl -= 360.

        if source.gl > self.gl_max or source.gl < self.gl_min:
            return False

        abs_gb = math.fabs(source.gb)
        if abs_gb > self.gb_max or abs_gb < self.gb_min:
            return False

        ra, dec = go.lb_to_radec(source.gl, source.gb)

        if ra > self.ra_max or ra < self.ra_min:
            return False
        if dec > self.dec_max or dec < self.dec_min:
            return False

        # Randomly decide if pulsar is in completed area of survey
        if random.random() > self.coverage:
            return False

        return True

    def dm_smear(self, source):
        """
        Calculate delay in pulse across a channel due to dm smearing. Formula's
        based on 'Handbook of Pulsar Astronomy" by Duncan Lorimer & Michael
        Kramer, section A2.4. Note the power of the forefactor has changed due
        to the central frequency being given in MHz.

        Args:
            source (class): Source object with a dm attribute
        Returns:
            t_dm, t_dm_err (float): Time of delay [ms] at central band
                                     frequency, with its error assuming a
                                     20% uncertainty in the dispersion measure
        """
        t_dm = 8.297616e6 * self.bw_chan * source.dm * (self.central_freq)**-3
        t_dm_err = t_dm / source.dm/(0.20*t_dm)
        return t_dm, t_dm_err

    def cal_flux(self, source, freq):
        """
        Calculate the flux of an FRB source at a particular frequency

        Args:
            source (class): Source method, needed for flux at 1400 MHz
            freq (float): Frequency [MHz] at which to calculate the flux
        Returns:
            flux (float): Source flux [TODO] at given input frequency
        """
        return source.s_1400() * (self.central_freq/freq)**source.spindex

    def cal_T_sky(self, source):
        """
        Calculate the sky temperature from the Haslam table, before scaling to
        the survey frequency. The temperature sky map is given in the weird
        units of HealPix and despite looking up info on this coordinate system,
        I don't have the foggiest idea of how to transform these to galactic
        coordinates. I have therefore directly copied the following code from
        psrpoppy in the assumption Sam Bates managed to figure it out.

        Args:
            source (class): Needed for coordinates
        Returns:
            T_sky (float): Sky temperature [K]
        """

        # ensure l is in range 0 -> 360
        B = source.gb
        if source.gl < 0.:
            L = 360 + source.gl
        else:
            L = source.gl

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

        return T_sky

    def calc_snr(self, source, population):
        """
        Calculate the signal to noise ratio of a source in the survey

        Args:
            source (class): Source of which to calculate the signal to noise
            population (class): Population to which the source belongs

        Returns:
            snr (float): Signal to noise ratio based on the radiometer equation
                         for a single pulse. Will return -2.0 if source is not
                         in survey region
        """

        if not self.in_region(source):
            return -2.0

        if self.gain_pattern == 'gaussian':

            # Formula's based on 'Interferometry and Synthesis in Radio
            # Astronomy' by A. Richard Thompson, James. M. Moran and
            # George W. Swenson, JR. (Second edition), around p. 15

            # Angular variable on sky, defined so that at fwhm/2, the
            # intensity profile is exactly 0.5. It's multiplied by the sqrt of
            # a random number to ensure the distribution of variables on the
            # sky remains uniform, not that it increases towards the centre
            # (uniform random point within circle). You could see this as
            # the calculated offset from the centre of the beam.
            xi = self.fwhm * math.sqrt(random.random()) / 2.

            # Intensity profile
            alpha = 2*math.sqrt(math.log(2))
            int_pro = math.exp(-(alpha*xi/self.fwhm)**2)

        # Dispersion measure across single channel, with error
        t_dm, t_dm_err = self.dm_smear(source)

        # Intrinsic pulse width
        w_int = source.width

        # Calculate scattering
        t_scat = go.scatter_bhat(source.dm, freq=self.central_freq)

        # From Narayan (1987, DOI: 10.1086/165442)
        # Also Cordes & McLaughlin (2003, DOI: 10.1086/378231)
        # For details see p. 30 of Emily Petroff's thesis (2016), found here:
        # http://hdl.handle.net/1959.3/417307
        w_eff = math.sqrt(w_int**2 + t_dm**2 + t_dm_err**2 +
                          t_scat**2 + self.t_samp**2)

        # Calculate total temperature
        T_sky = self.cal_T_sky(source)
        T_tot = self.T_sys + T_sky

        # Calculate flux density at central frequency
        s_peak = source.s_1400()

        # Radiometer equation for single pulse (Dewey et al., 1984)
        snr = s_peak * self.gain * math.sqrt(self.n_pol*self.bw_chan*w_eff)
        snr /= (T_tot * self.beta)

        # Account for offset in beam
        snr *= int_pro

        return snr

    def scint(self, source, snr):
        """
        Calculate scintillation effect on the signal to noise ratio (rather than
        adapting the flux, as the snr can change per survey attempt). Formula's
        based on 'Handbook of Pulsar Astronomy" by Duncan Lorimer & Michael
        Kramer, section 4.2.

        Args:
            src (class): Source object
            snr (float): Signal to noise ratio
        Returns:
            snr (float): Signal to noise ratio modulated by scintillation
        """
        # Calculate scattering
        t_scat = go.scatter_bhat(source.dm, freq=self.central_freq)
        # Convert to seconds
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

            t_diss, decorr_bw = go.ne2001_scint_time_bw(source.dist,
                                                        source.gl,
                                                        source.gb,
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
        snr = random.gauss(snr, m*snr)

        return snr

