import os

from collections import OrderedDict as OD

from log import Log


class Population:
    """Class to hold a population of FRBs"""

    def __init__(self,
                 electron_model=None,
                 lum_dist_pars=None,
                 log_loc=None,
                 no_log=False,
                 quiet=False,
                 verbose=False):

        # Population properties
        self.electron_model = None
        self.lum_dist_pars = None

        # Store FRB sources
        self.sources = []

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

    def write_csv(self, out=None, sep=','):
        """
        Write out source properties as csv file

        Args:
            out (str): Outfile location (default=data/logs/population.csv)
            sep (str): [Optional] Define seperator
        """

        if out is None:
            loc = '../data/logs/population.csv'
            out = os.path.join(os.path.dirname(__file__), loc)

        with open(out, 'w') as f:

            # Find all source properties
            a = self.sources[0].__dict__
            attrs = OD(sorted(a.items()))

            # Create header
            text = sep.join(attrs.keys()) + '\n'

            # Print values per source
            for src in self.sources:
                text += sep.join([str(src.__dict__[k]) for k in attrs]) + '\n'

            f.write(text)

    def write_ascii(self, out=None):
        """
        Write out source properties in ascii format

        Args:
            out (str): Outfile location (default=data/logs/population.dat)
        """

        if out is None:
            loc = '../data/logs/population.dat'
            out = os.path.join(os.path.dirname(__file__), loc)

        self.write_csv(out=out, sep=' ')
