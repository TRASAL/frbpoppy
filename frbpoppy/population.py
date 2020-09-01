"""Define a class to hold a population of FRBs."""
import os
import dill as pickle
import numpy as np
from copy import deepcopy

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
        output = open(path, 'wb')
        pickle.dump(self, output, 2)
        output.close()

    def n_sources(self):
        """Return the number of FRB sources."""
        return len(self.frbs.ra)

    def n_srcs(self):
        return self.n_sources()

    def n_bursts(self):
        """Return the number of bursts."""
        try:  # Will only work for a repeater population
            n = np.count_nonzero(~np.isnan(self.frbs.time))
        except TypeError:
            n = self.n_sources()
        return n

    def n_repeaters(self):
        """Return the numer of repeaters in a population."""
        try:
            return np.sum((~np.isnan(self.frbs.time)).sum(1) > 1)
        except TypeError:
            return 0

    def n_rep(self):
        return self.n_repeaters()

    def n_one_offs(self):
        """Return the numer of one-offs in a population."""
        try:
            return np.sum((~np.isnan(self.frbs.time)).sum(1) <= 1)
        except TypeError:
            return self.n_sources()

    def n_oneoffs(self):
        """Return the numer of one-offs in a population."""
        return self.n_one_offs()


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


def split_pop(pop, mask=None):
    """Split a population.

    Args:
        pop (Population): Population to be split
        mask (array): Numpy boolean mask. If none, split into repeaters and
            one_offs in that order

    Returns:
        tuple: Tuple of population classes

    """
    if mask is None:
        mask = ((~np.isnan(pop.frbs.time)).sum(1) > 1)
    pop_true = deepcopy(pop)
    pop_false = deepcopy(pop)
    pop_true.frbs.apply(mask)
    pop_false.frbs.apply(~mask)
    return pop_true, pop_false


def merge_pop(*args, random=False):
    """Merge populations.

    Args:
        Populations to merge
        random (bool): If wishing to shuffle the frbs from different pops

    Returns:
        Population

    """
    mp = args[0]  # Main population

    for pop in args:
        # Merge each parameter
        for attr in mp.frbs.__dict__.keys():
            parm = getattr(mp.frbs, attr)
            if type(parm) is np.ndarray:
                parms = []
                for pop in args:
                    parms.append(getattr(pop.frbs, attr))

                try:
                    merged_parm = np.concatenate(parms, axis=0)
                except ValueError:
                    # Check maximum size values should be padded to
                    max_size = max([p.shape[1] for p in parms])
                    new_parms = []

                    # Ensure matrices are the same shapes by padding them
                    for p in parms:
                        if p.shape[1] != max_size:
                            padded_p = np.zeros((p.shape[0], max_size))
                            padded_p[:] = np.nan
                            padded_p[:, :p.shape[1]] = p
                            new_parms.append(padded_p)
                        else:
                            new_parms.append(p)

                    merged_parm = np.concatenate(new_parms, axis=0)

                setattr(mp.frbs, attr, merged_parm)

    if random:
        shuffle = np.random.permutation(mp.frbs.z.shape[0])
        for attr in mp.frbs.__dict__.keys():
            parm = getattr(mp.frbs, attr)
            if type(parm) is np.ndarray:
                setattr(mp.frbs, attr, parm[shuffle])

    return mp
