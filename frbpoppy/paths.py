"""Set all possible needed paths."""
import os.path
from frbpoppy.log import pprint


class Paths():
    """All directory paths."""

    def __init__(self):
        """Initialise."""
        # Convenient paths to have
        self.code = os.path.dirname(__file__)
        self.home = os.path.expanduser("~")
        self.downloads = os.path.expanduser("~/Downloads")

        self.subfolders = ['data',
                           'frbcat',
                           'results',
                           'populations',
                           'surveys',
                           'models']

        self.config = {s: '' for s in self.subfolders}

    def __str__(self):
        """Define how to print a path object to a console."""
        s = 'Paths:'

        attributes = []
        for e in self.subfolders:
            attr = '\n\t{0:12.11}{1:.60}'.format(e, getattr(self, e)())
            attributes.append(attr)

        s += ''.join(attributes)

        return s

    def check(self, path):
        """Perform checks on path."""
        # Just convient to have files ending in a slash
        if path[-1] != '/':
            path += '/'

        if not os.path.exists(path):
            pprint(f"Creating directory {path}")
            os.makedirs(path)

        return path

    def store(self, name, default, *args):
        """
        Where files are to be stored.

        Args:
            name (str): Which group of files, populations, results etc
            default (str): Default path for the group

        Returns:
            path (str): The path of the group

        """
        # Check whether a new default is being set
        if args:
            self.config[name] = args[0]

        # If a new default has been set
        if self.config[name]:
            path = self.config[name]
        else:
            path = default
        # Create directory if not present
        path = self.check(path)

        return path

    def data(self, *args):
        """Where all results are to be stored."""
        default = os.path.realpath(os.path.join(self.code, '../data/')) + '/'
        return self.store('data', default, *args)

    def results(self, *args):
        """Where all results are to be stored."""
        default = self.data() + 'results/'
        return self.store('results', default, *args)

    def populations(self, *args):
        """Where all populations are to be stored."""
        default = self.results()
        return self.store('populations', default, *args)

    def surveys(self, *args):
        """Where all surveys are to be stored."""
        default = self.data() + 'surveys/'
        return self.store('surveys', default, *args)

    def models(self, *args):
        """Where all models are to be stored."""
        default = self.data() + 'models/'
        return self.store('models', default, *args)

    def frbcat(self, *args):
        """Where all frbcatalogues are to be stored."""
        default = self.data() + 'frbcat/'
        return self.store('frbcat', default, *args)


paths = Paths()
