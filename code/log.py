import inspect
import logging
import logging.config
import os
import sys

def pprint(*s):
    """Hack to make for more informative print statements"""
    f = inspect.stack()[1][1].split('/')[-1]
    m = '{:13.13} |'.format(f)
    print(m, *s)


class Log(object):
    """Setting default logging styles. Only really helpful when invoking as

    >>> logger = Log().logger()

    and to log/output text, you can use

    >>> logger.info('Something moderately important')

    Args:
        save (boolean): Log to file
        verbose (boolean): Give more output in console
        quiet (boolean): Stop console output
        loc (str): Give the log's filename
    """

    def __init__(self,
                 no_log=False,
                 verbose=False,
                 quiet=False,
                 loc=None):
        # Input
        self.save = not no_log
        self.verbose = verbose
        self.quiet = quiet
        self.loc = loc

        # Config options
        self.handlers = {}
        self.formatters = {}
        self.level = logging.INFO
        # Find name of the executed program
        self.filename = sys.argv[0].split('/')[-1]

        # Apply input
        self._file()
        self._console()
        self._cfg()

    def _file(self):
        """Add file options"""

        # To save or not to save
        if not self.save:
            return

        # Set default location of log file
        if self.loc is None:
            fn = self.filename.split('.')[-2]
            self.loc = __file__[:-6] + '../data/logs/' + fn + '.log'

        # Warn if log size is getting out of hand
        try:
            if os.path.getsize(self.loc) > 1e6:
                print('WARNING: Your log file is getting rather big - ' +
                      'please consider deleting it')
        except FileNotFoundError:
            pass

        # Format log file output
        frmt_file = '%(asctime)s | %(name)-12s | %(levelname)-8s | %(message)s'
        self.formatters['f'] = {'format': frmt_file}

        hnd_file = {'class': 'logging.handlers.RotatingFileHandler',
                    'formatter': 'f',
                    'level': logging.DEBUG,
                    'filename': self.loc,
                    'mode': 'a',
                    'backupCount': 0}  # Increase for longer backups
        self.handlers['file'] = hnd_file

    def _console(self):
        """Add console options"""

        if self.verbose and self.quiet:
            raise ValueError("Really now.. You can't ask me to be " +
                             "both verbose and quiet")

        frmt_style = 'c'

        if self.verbose:
            self.level = logging.DEBUG
            frmt_style = 'f'
        if self.quiet:
            self.level = logging.ERROR

        # Format console output
        frmt_console = '%(filename)-12s | %(message)s'
        self.formatters['c'] = {'format': frmt_console}

        hnd_console = {'class': 'logging.StreamHandler',
                       'formatter': frmt_style,
                       'level': self.level}
        self.handlers['console'] = hnd_console

    def _cfg(self):
        """Setup configuration options"""

        cfg = dict(version=1,
                   disable_existing_loggers=False,
                   formatters=self.formatters,
                   handlers=self.handlers,
                   root={'handlers': [k for k in self.handlers],
                         'level': logging.DEBUG}
                   )

        logging.config.dictConfig(cfg)

    def logger(self):
        """Return logger instance (see argparse library)"""

        # Finds out which function/class is calling this log class
        module = inspect.stack()[1][1].split('/')[-1]
        if module == '<module>':
            module = self.filename

        return logging.getLogger(module)
