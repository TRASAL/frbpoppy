from collections import OrderedDict as OD
from io import StringIO
import os
import pickle
import pandas as pd

from log import pprint


class Population:
    """Class to hold a population of FRBs"""

    def __init__(self):

        # Population properties
        self.cosmology = None
        self.electron_model = None
        self.f_max = None
        self.f_min = None
        self.H_0 = None
        self.lum_max = None
        self.lum_min = None
        self.lum_pow = None
        self.name = None
        self.repeat = None
        self.si_mean = None
        self.si_sigma = None
        self.time = None
        self.v_max = None
        self.w_max = None
        self.w_min = None
        self.W_m = None
        self.W_v = None
        self.z_max = None

        # Store FRB sources
        self.sources = []
        self.n_srcs = 0

    def __str__(self):
        """Define how to print a population object to a console."""

        s = 'Population properties:'

        attributes = []
        for e in self.__dict__:
            attr = '\n\t{0:12.11}{1:.60}'.format(e, str(self.__dict__[e]))
            attributes.append(attr)

        s += ''.join(attributes)

        return s

    def add(self, source):
        """Add a source to the population."""
        self.sources.append(source)
        self.n_srcs += 1

    def save(self, out=None, sep=','):
        """
        Write out source properties as data file.

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
        src_attrs = sorted(list(a.keys()), key=lambda v: v.upper())
        src_attrs.remove('detected')
        src_attrs.remove('frbs')

        # Add frb properties
        a = self.sources[0].frbs[0].__dict__
        frb_attrs = sorted(list(a.keys()), key=lambda v: v.upper())
        frb_attrs.remove('detected')

        # Create header
        attrs = src_attrs + frb_attrs
        data = sep.join(attrs) + '\n'

        # Print values per source
        for src in self.sources:

            src_data = [str(src.__dict__[k]) for k in src_attrs]

            for frb in src.frbs:
                frb_data = [str(frb.__dict__[k]) for k in frb_attrs]
                data += sep.join(src_data + frb_data) + '\n'

        return data

    def to_df(self):
        """Gather source values into a pandas dataframe"""
        data = StringIO(self.values())
        df = pd.read_csv(data)
        return df

    def pickle_pop(self, out=None):
        """Allow the population to be pickled for future use."""
        # Check if an output file has been given
        if out is None:
            loc = '../data/results/population_{}.p'.format(self.name.lower())
            out = os.path.join(os.path.dirname(__file__), loc)

        output = open(out, 'wb')
        pickle.dump(self, output, 2)
        output.close()


def unpickle(filename=None):
    """Quick function to unpickle a population

    Args:
        filename (str, optional): Define the path to the pickled population,
            or give the population name

    Returns:
        pop (Population): Population class
    """

    # Find population file
    if os.path.isfile(filename):
        f = open(filename, 'rb')
    else:
        # Find standard population files
        try:
            p = '../data/results/population_{}.p'.format(filename.lower())
            path = os.path.join(os.path.dirname(__file__), p)
            f = open(path, 'rb')
        except FileNotFoundError:
            s = 'Pickled population file "{0}" does not exist'.format(filename)
            raise FileNotFoundError(s)

    # Unpickle
    pop = pickle.load(f)
    f.close()

    return pop
