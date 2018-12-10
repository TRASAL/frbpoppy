"""Classes to hold rate counters."""
from copy import deepcopy
from frbpoppy import pprint


class Rates:
    """Class to hold rate counters."""

    def __init__(self):
        """Initializing."""
        # Rates
        self.det = 0  # Number detected
        self.faint = 0  # Number too faint to detect
        self.late = 0  # Number too late to detect
        self.out = 0  # Number outside survey, space or timewise
        self.vol = 0  # Number per Gpc^3
        self.days = 0  # Days of a survey
        self.name = ''  # Name of a survey

        # Scaling factors
        self.f_area = 1  # Scaling factor area
        self.f_time = 1  # Scaling factor time
        self.scaled_area = False
        self.scaled_time = False

    def __str__(self):
        """How to print the class."""
        # Set up title
        r = '{:20.19} {:>10} {:>10}\n'
        t = r.format(self.name, 'Days', 'FRBs')
        line = '-'*len(t.split('\n')[-2].strip()) + '\n'
        t += line

        # Format rates
        rdays = round(self.days)
        t += r.format('In population', rdays, round(self.tot()))
        t += r.format('Detected', rdays, round(self.det))
        t += r.format('Too late', rdays, round(self.late))
        t += r.format('Too faint', rdays, round(self.faint))
        t += r.format('Outside survey', rdays, round(self.out))
        t += r.format('/Gpc^3', 365.25, round(self.vol))
        t += r.format('Expected', round(self.exp, 4), 1)
        t += line

        return pprint(t, output=False)

    def tot(self):
        """Calculate the total number of rates."""
        return self.det + self.out + self.faint + self.late

    @property
    def exp(self):
        """Days before an FRB is detected."""
        try:
            return self.days/self.det
        except ZeroDivisionError:
            return float('NaN')


def scale(rates, area=True, time=False):
    """Scale rates."""
    rates = deepcopy(rates)

    if area:
        rates.scaled_area = True
        tot = rates.tot()
        rates.det *= rates.f_area
        rates.late *= rates.f_area
        rates.faint *= rates.f_area
        rates.out = tot - rates.det - rates.faint - rates.late

    if time:
        rates.scaled_time = True
        rates.det *= rates.f_time
        rates.late *= rates.f_time
        rates.faint *= rates.f_time
        rates.out *= rates.f_time
        rates.days = 1

    return rates
