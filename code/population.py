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
        self.n_srcs = 0
        self.name = None  # Population name

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

    def save(self, out=None, sep=','):
        """
        Write out source properties as data file

        Args:
            out (str): Outfile location (default=data/results/population.csv or
                       data/results/population_<survey_name>.csv if survey)
            sep (str): Define seperator in file, which also changes the file
                       type between .dat and .csv
        """
        # Check the population contains sources
        if len(self.sources) == 0:
            print('Nothing to write: survey population contains no sources')
            return

        # Set default file locations
        if out is None:

            # Check if a population has been a survey name
            if self.name is None:
                loc = '../data/results/population'
                out = os.path.join(os.path.dirname(__file__), loc)
            else:
                loc = '../data/results/population_' + self.name.lower()
                out = os.path.join(os.path.dirname(__file__), loc)

            # Set file types
            if sep == ',':
                out += '.csv'
            elif sep == ' ':
                out += '.dat'

        with open(out, 'w') as f:
            f.write(self.values(sep=sep))

    def values(self, sep=','):
        """Gather source values into table"""

        # Check the population contains sources
        if len(self.sources) == 0:
            print('Nothing to write: survey population contains no sources')
            return

        # Find all source properties
        a = self.sources[0].__dict__
        attrs = OD(sorted(a.items()))

        # Create header
        data = sep.join(attrs.keys()) + '\n'

        # Print values per source
        for src in self.sources:
            data += sep.join([str(src.__dict__[k]) for k in attrs]) + '\n'
