"""Define a class to hold a population of FRBs."""
from io import StringIO
import os
import pickle
import pandas as pd

from frbpoppy.log import pprint
from frbpoppy.paths import paths


class Population:
    """Class to hold a population of FRBs."""

    def __init__(self):
        """Initializing."""
        # Population properties
        self.name = None

        # Frequency emission limits [MHz]
        self.f_max = None
        self.f_min = None

        # Dispersion Measure [pc/cm^3]
        self.dm_host_model = None
        self.dm_host_mu = None
        self.dm_host_sigma = None
        self.dm_igm_index = None
        self.dm_igm_sigma = None
        self.dm_mw_model = None

        # Luminosity
        self.lum_max = None
        self.lum_min = None
        self.lum_pow = None

        # Spectral index
        self.si_mu = None
        self.si_sigma = None

        # Pulse width
        self.w_model = None
        self.w_mu = None
        self.w_sigma = None
        self.w_max = None
        self.w_min = None

        # Cosmology
        self.dist_co_max = None
        self.H_0 = None
        self.n_model = None
        self.vol_co_max = None
        self.W_m = None
        self.W_v = None
        self.z_max = None

        # Repeater properties
        self.repeat = None
        self.time = None  # seconds

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

    def get(self, par):
        """Get a list of a parameter."""
        pars = []

        for src in self.sources:
            try:
                pars.append(getattr(src, par))
            except AttributeError:
                for frb in src.frbs:
                    pars.append(getattr(frb, par))

        return pars

    def to_df(self):
        """Gather source values into a pandas dataframe."""
        values = self._values()

        if not values:
            return None

        if len(values) >= 200000:
            # Take quarter of the values
            pprint(f'Quartering {self.name} population')
            values = values[:len(values)//4]
            index = values.rfind('\n')
            values = values[:index]

        data = StringIO(values)
        df = pd.read_csv(data)
        return df

    def save(self, extention='p'):
        """
        Write out source properties as data file.

        Args:
            extention (str): Type of file to save.
                Choice from ['csv', 'dat', 'p']
        """
        # Check if a population has been a survey name
        if not self.name:
            file_name = 'population'
        else:
            file_name = 'population_' + self.name.lower()

        path = paths.populations() + file_name

        # Set file types
        if extention == 'p':
            path += '.p'
            self.to_pickle(path)
        if extention == 'csv':
            path += '.csv'
            self._to_csv(path)
        elif extention == 'dat':
            path += '.dat'
            self._to_csv(path, sep=' ')

    def _to_csv(self, path, sep=','):
        """Write a population to a csv file.

        Args:
            path (str): Path to which to write
            sep (str): Seperator character

        """
        with open(path, 'w') as f:
            v = self._values(sep=sep)
            if not v:
                v = ' '
            f.write(v)

    def _values(self, sep=','):
        """
        Gather source values into table in string format.

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

    def to_pickle(self, path):
        """Write a population to a pickled file for future use.

        Args:
            path (str): Path to which to write

        """
        output = open(path, 'wb')
        pickle.dump(self, output, 2)
        output.close()


def unpickle(filename=None):
    """Quick function to unpickle a population.

    Args:
        filename (str, optional): Define the path to the pickled population,
            or give the population name

    Returns:
        Population: Population class

    """
    # Find population file
    if os.path.isfile(filename):
        f = open(filename, 'rb')
    else:
        # Find standard population files
        try:
            name = filename.lower()
            p = paths.populations() + f'population_{name}.p'
            f = open(p, 'rb')
        except FileNotFoundError:
            s = 'Pickled population file "{0}" does not exist'.format(filename)
            raise FileNotFoundError(s)

    # Unpickle
    pop = pickle.load(f)
    f.close()

    return pop
