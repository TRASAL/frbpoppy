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
        self.gain_pattern = pattern
        self.aa = False  # Whether aperture array

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
        Attempt to parse survey files

        Args:
            f (str): Filename, see Survey class
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
                self.t_sys = float(v)  # [K]
            elif p.count('centre frequency'):
                self.freq = float(v)  # [MHz]
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

    def calc_snr(self, source, population):
        """
        Calculate the signal to noise ratio of a source in the survey

        Args:
            source (class): Source of which to calculate the signal to noise
            population (class): Population to which the source belongs

        Returns:
            0.0: If source is not in survey region
        """

        if not self.in_region(source):
            return 0.0

        if self.gain_pattern == 'gaussian':

            # Formula's based on 'Interferometry and Synthesis in Radio
            # Astronomy' by A. Richard Thompson, James. M. Moran and
            # George W. Swenson, JR. (Second edition), around p. 15

            # Angular variable on sky, defined so that at fwhm/2, the
            # intensity profile is exactly 1. It's multiplied by the sqrt of
            # a random number to ensure the distribution of variables on the
            # sky remains uniform, not that it increases towards the centre
            # (uniform random point within circle)
            xi = self.fwhm * math.sqrt(random.random()) / 2.

            # Intensity profile
            alpha = 2*math.sqrt(math.log(2))
            int_pro = math.exp(-(alpha*xi/self.fwhm)**2)

