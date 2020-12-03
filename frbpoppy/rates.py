"""Classes to hold rate counters."""
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
        f = '{:20.19} {:>10} {:>10} {:>10}\n'
        t = f.format(self.name, 'Days', f'{self.object_type.title()}s', '%')
        line = '-'*len(t.split('\n')[-2].strip()) + '\n'
        t += line

        def r(value, d=4):
            """Round a value"""
            return round(value, d)

        def per(value):
            """Calculate the percentage."""
            return r(value/self.tot * 100)

        # Format rates
        days = r(self.days)
        t += f.format('Cosmic Population', days, r(self.tot), 100)
        t += f.format('Too late', days, r(self.late), per(self.late))
        t += f.format('Outside survey', days, r(self.out), per(self.out))
        t += f.format('Outside pointings', days, r(self.pointing),
                      per(self.pointing))
        t += f.format('Too faint', days, r(self.faint), per(self.faint))
        t += f.format('Detected', days, r(self.det), per(self.det))
        t += f.format('/Gpc^3', 365.25, r(self.vol, 2), '-')
        t += f.format('Expected', r(self.exp, 4), 1, '-')
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
