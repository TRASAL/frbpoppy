"""Set all possible needed paths."""
import os.path


class FileObject():
    """Either file or directory.

    Examples:
        >>> pop_path = FileObject()
        >>> pop_path.set_dir('/random/path/')
        /random/path/
        >>> print(pop_path.dir)
        /random/path/
        >>> pop_path.set_file = 'population_test.csv'
        /random/path/population_test.csv
        >>> print(pop_path.file)
        /random/path/population_test.csv

    """

    def __init__(self):
        """Initialise."""
        self.dir = ''  # Directory path
        self.file = ''  # Path to file

    def set_dir(self, dir_name):
        """Check whether directory exsists, if not, will create it."""
        self.dir = dir_name
        if not os.path.exists(self.dir):
            print(f"Files will be saved in {self.dir}")
            os.makedirs(self.dir)
        return self.dir

    def set_file(self, file_name):
        """Return full path to file with file_name."""
        self.file = os.path.join(self.dir, file_name)
        return self.file


class Path():
    """All directory paths."""

    def __init__(self):
        """Initialise."""
        # Convenient paths to have
        self._here = os.path.dirname(__file__)
        self._home = os.path.expanduser("~")

        # Set up where to save the results
        self.downloads = os.path.expanduser("~/Downloads")
        self.output = FileObject()
        self.output.set_dir(os.path.join(self.downloads, 'frbpoppy/'))

        self.results = FileObject()
        self.results.set_dir(self.rel('results'))

        self.frbcat = FileObject()
        self.frbcat.set_dir(self.rel('frbcat'))

        self.logs = FileObject()
        self.logs.set_dir(self.rel('logs'))

        self.models = FileObject()
        self.models.set_dir(self.rel('models'))

        self.surveys = FileObject()
        self.surveys.set_dir(self.rel('surveys'))

    def rel(self, path, to=None):
        """Return relative path."""
        if not to:
            to = self.output.dir
        return os.path.join(to, path) + '/'


path = Path()
