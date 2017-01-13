import sys

from collections import OrderedDict as OD

from log import Log
from source import Source

class Population:
    """Class to hold a population of FRBs"""


    def __init__(self,
                 log_loc=None,
                 no_log=False,
                 quiet=False,
                 verbose=False):

        # Store FRB source
        self.population = []

        # Logging options
        self.log_loc = log_loc
        self.no_log = no_log
        self.quiet = quiet
        self.verbose = verbose
        self.logger = self.log()

        # Counter
        self.n_det = 0


    def __str__(self):
        """Define how to print a population object to a console"""

        s = 'Population properties:'

        attributes = []
        for e in self.__dict__:
            attr = '\n\t{0:11.12}{1:.60}'.format(e, str(self.__dict__[e]))
            attributes.append(attr)

        s += ''.join(attributes)

        return s


    def log(self):
        """Set up log"""
        logger = Log(no_log=self.no_log,
                     verbose=self.verbose,
                     quiet=self.quiet,
                     loc=self.log_loc).logger()
        return logger


    def write_csv(self,out):
        """Write out source properties as csv file"""

        with open(out, 'w') as f:

            # Find all source properties
            a = self.population[0].__dict__
            attrs = OD(sorted(a.items()))

            # Create header
            text = '#' + ','.join(attrs.keys()) + '\n'

            # Print values per source
            for src in self.population:
                text += ','.join([str(src.__dict__[k]) for k in attrs]) + '\n'

            f.write(text)


    def write_pickle(self,out):
        """Write a pickled class"""
        return 'Hmm'

