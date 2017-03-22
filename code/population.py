import os
from collections import OrderedDict as OD

from log import pprint


class Population:
    """Class to hold a population of FRBs"""

    def __init__(self,
                 cosmology=None,
                 electron_model=None,
                 f_max=None,
                 f_min=None,
                 H_0=None,
                 lum_max=None,
                 lum_min=None,
                 lum_pow=None,
                 name=None,
                 si_mean=None,
                 si_sigma=None,
                 v_max=None,
                 W_m=None,
                 W_v=None,
                 z_max=None):

        # Population properties
        self.cosmology = cosmology
        self.electron_model = electron_model
        self.f_max = f_max
        self.f_min = f_min
        self.H_0 = H_0
        self.lum_max = lum_max
        self.lum_min = lum_min
        self.lum_pow = lum_pow
        self.name = name
        self.si_mean = si_mean
        self.si_sigma = si_sigma
        self.v_max = v_max
        self.W_m = W_m
        self.W_v = W_v
        self.z_max = z_max

        # Store FRB sources
        self.sources = []

        # Counter
        self.n_srcs = 0

    def __str__(self):
        """Define how to print a population object to a console"""

        s = 'Population properties:'

        attributes = []
        for e in self.__dict__:
            attr = '\n\t{0:12.11}{1:.60}'.format(e, str(self.__dict__[e]))
            attributes.append(attr)

        s += ''.join(attributes)

        return s

    def save(self, out=None, sep=','):
        """Write out source properties as data file

        Args:
            out (str): Outfile location. Defaults to
                data/results/population.csv or
                data/results/population_<survey_name>.csv if survey
            sep (str): Define seperator in file, which also changes the file
                       type between .dat and .csv. Defaults to csv
        """

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
            v = self.values(sep=sep)
            if not v:
                v = ' '
            f.write(v)

    def values(self, sep=','):
        """Gather source values into table in string format

        Args:
            sep (str, optional): Define the seperator between values

        Returns:
            data (str): Data table with the values of all attributes, including
                a header at the top
        """

        # Check the population contains sources
        if len(self.sources) == 0:
            m = 'Population {} contains no sources'
            pprint(m.format(self.name))
            return

        # Find all source properties
        a = self.sources[0].__dict__
        attrs = OD(sorted(a.items()))

        # Create header
        data = sep.join(attrs.keys()) + '\n'

        # Print values per source
        for src in self.sources:
            data += sep.join([str(src.__dict__[k]) for k in attrs]) + '\n'

        return data
