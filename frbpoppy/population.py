"""Define a class to hold a population of FRBs."""
import os
import dill as pickle
import numpy as np
from copy import deepcopy
from types import LambdaType

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

        # Store FRB sources
        self.frbs = FRBs()
        self.uid = None  # Unique Identifier

    def __str__(self):
        """Define how to print a population object to a console."""
        s = 'Population properties:'

        # TODO: Code this to print all properties

        return s

    def to_df(self):
        """Gather source values into a pandas dataframe."""
        df = self.frbs.to_df()
        return df

    def save(self, path=None):
        """
        Write out source properties as data file.

        Args:
            path (str): Path to which to save.
        """
        if path is None:
            # Check if a population has been a survey name
            if not self.name:
                file_name = 'pop'
            else:
                file_name = self.name.lower()

            if self.uid:
                file_name += f'_{self.uid}'

            path = paths.populations() + f'{file_name}.p'

        self.to_pickle(path)

    def to_csv(self, path):
        """Write a population to a csv file.

        Args:
            path (str): Path to which to write

        """
        df = self.frbs.to_df()
        df.to_csv(path)

    def to_pickle(self, path):
        """Write a population to a pickled file for future use.

        Args:
            path (str): Path to which to write

        """
        # # Python doesn't support the pickling of lambda functions
        # def islambda(obj):
        #     """Check whether variable is lambda instance."""
        #     obj_type = isinstance(obj, LambdaType)
        #     return obj_type and obj.__name__ == '<lambda>'
        #
        # # Rename lambda functions
        # for attr in self.__dict__.keys():
        #     pprint(attr)
        #     if attr.endswith('func'):
        #         import IPython; IPython.embed()
        #     parm = getattr(self, attr)
        #     if islambda(parm):
        #         parm.__name__ = attr
        #         setattr(self, attr, parm)

        output = open(path, 'wb')
        pickle.dump(self, output, 2)
        output.close()

    def n_sources(self):
        """Return the number of FRB sources."""
        return len(self.frbs.ra)

    def n_bursts(self):
        """Return the number of bursts."""
        try:  # Will only work for a repeater population
            n = np.count_nonzero(~np.isnan(self.frbs.time))
        except TypeError:
            n = self.n_sources()
        return n


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
