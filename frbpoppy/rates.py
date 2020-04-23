"""Classes to hold rate counters."""
from copy import deepcopy
import numpy as np
from frbpoppy.misc import pprint


class Rates:
    """Class to hold rate counters."""

    def __init__(self, object_type='source'):
        """Initialize.

        Args:
         object_type (str): Type of object of which to keep track.
        """
        # Rates
        self.object_type = object_type
        self.det = 0  # Number detected
        self.faint = 0  # Number too faint to detect
        self.late = 0  # Number too late to detect
        self.out = 0  # Number outside survey region
        self.pointing = 0  # Number outside of pointing
        self.vol = 0  # Number per Gpc^3
        self.days = 0  # Days of a survey
        self.name = ''  # Name of a survey
        self.tot = 0  # Total number of detectable objects

        # Scaling factors
        self.f_area = 1  # Scaling factor area
        self.scaled_area = False

    def __str__(self):
        """How to print the class."""
        # Set up title
        r = '{:20.19} {:>10} {:>10}\n'
        t = r.format(self.name, 'Days', f'{self.object_type.title()}s')
        line = '-'*len(t.split('\n')[-2].strip()) + '\n'
        t += line

        # Format rates
        rdays = round(self.days, 3)
        t += r.format('Cosmic Population', rdays, round(self.tot))
        t += r.format('Too late', rdays, round(self.late, 3))
        t += r.format('Outside survey', rdays, round(self.out, 3))
        t += r.format('Outside pointings', rdays, round(self.pointing, 3))
        t += r.format('Too faint', rdays, round(self.faint, 3))
        t += r.format('Detected', rdays, round(self.det, 3))
        t += r.format('/Gpc^3', 365.25, round(self.vol, 2))
        t += r.format('Expected', round(self.exp, 4), 1)
        t += line

        return pprint(t, output=False)

    def sum(self):
        """Calculate the total number of rates."""
        return self.out + self.late + self.pointing + self.faint + self.det

    @property
    def exp(self):
        """Days before an FRB is detected."""
        try:
            return self.days/self.det
        except ZeroDivisionError:
            return float('NaN')

    def scale_by_area(self):
        """Scale rates."""
        self.scaled_area = True
        tot = self.tot
        self.det *= self.f_area
        self.late *= self.f_area
        self.faint *= self.f_area
        self.out = np.abs(tot - self.det - self.faint - self.late)
