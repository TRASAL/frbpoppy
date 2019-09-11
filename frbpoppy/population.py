"""Define a class to hold a population of FRBs."""
from io import StringIO
import os
import pickle
import numpy as np
import pandas as pd
from copy import deepcopy

from frbpoppy.log import pprint
from frbpoppy.paths import paths
from frbpoppy.frbs import FRBs


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

        # Pulse width [ms]
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

        # Store FRB sources
        self.frbs = FRBs()
        self.n_frbs = 0
        self.uid = None  # Unique Identifier

    def __str__(self):
        """Define how to print a population object to a console."""
        s = 'Population properties:'

        # TODO: Code this to print all properties

        return s

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
            file_name = 'pop'
        else:
            file_name = self.name.lower()

        if self.uid:
            file_name += f'_{self.uid}'

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

    def to_csv(self, path, sep=','):
        """Write a population to a csv file.

        Args:
            path (str): Path to which to write
            sep (str): Seperator character

        """
        df = self.frbs.to_df(sep=sep)
        df.to_csv(path)

    def to_pickle(self, path):
        """Write a population to a pickled file for future use.

        Args:
            path (str): Path to which to write

        """
        output = open(path, 'wb')
        pickle.dump(self, output, 2)
        output.close()


def unpickle(filename=None, uid=None):
    """Quick function to unpickle a population.

    Args:
        filename (str, optional): Define the path to the pickled population,
            or give the population name
        uid (str, optional): Unique Identifier

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
            if uid:
                name += f'_{uid}'
            p = paths.populations() + f'{name}.p'
            f = open(p, 'rb')
        except FileNotFoundError:
            s = 'Pickled population file "{0}" does not exist'.format(filename)
            raise FileNotFoundError(s)

    # Unpickle
    pop = pickle.load(f)
    f.close()
    return pop


def split_pop(pop, mask):
    """Split a population.

    Args:
        pop (Population): Population to be split
        mask (Numpy): Numpy boolean mask

    Returns:
        tuple: Tuple of population classes

    """
    pop_true = deepcopy(pop)
    pop_false = deepcopy(pop)
    pop_true.frbs.apply(mask)
    pop_false.frbs.apply(~mask)
    return pop_true, pop_false
