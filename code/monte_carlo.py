"""Monte Carlo over various frbpoppy parameters."""

from collections import namedtuple
import numpy as np
import os
import pandas as pd
import sqlite3


class MonteCarlo:

    def __init__(self):

        self.Par = namedtuple('par', 'min max step set')

        self.dm_host = self.Par(0, 200, 10, 100)
        self.dm_igm_slope = self.Par(1000, 1400, 20, 1200)
        self.freq_min = self.Par(10e5, 10e10, 0.5e1, 10e6)
        self.freq_max = self.Par(10e5, 10e10, 0.5e1, 10e9)
        self.lum_bol_slope = self.Par(0.5, 1.5, 0.1, 1.)
        self.lum_bol_min = self.Par(1e30, 1e60, 1e1, 1e40)
        self.lum_bol_max = self.Par(1e30, 1e60, 1e1, 1e50)
        self.n_day = self.Par(2000, 14000, 2000, 10000)
        self.rep = self.Par(0.0, 0.1, 0.01, 0.05)
        self.si_mean = self.Par(-2.0, -1, 0.1, -1.4)
        self.si_sigma = self.Par(0.0, 0.5, 0.1, 0.0)
        self.w_int_min = self.Par(0.0, 5, 1.0, 1.)
        self.w_int_max = self.Par(0.0, 5, 1.0, 5.)

        self.pars = {'dm_host': self.dm_host,
                     'dm_igm_slope': self.dm_igm_slope,
                     'freq_min': self.freq_min,
                     'freq_max': self.freq_max,
                     'lum_bol_slope': self.lum_bol_slope,
                     'lum_bol_min': self.lum_bol_min,
                     'lum_bol_max': self.lum_bol_max,
                     'n_day': self.n_day,
                     'rep': self.rep,
                     'si_mean': self.si_mean,
                     'si_sigma': self.si_sigma,
                     'w_int_min': self.w_int_min,
                     'w_int_max': self.w_int_max}

        self.dfs = {}

    def par_range(self, p, mi=None, ma=None, st=None, se=None):
        """Quick range generator."""
        if not mi:
            mi = p.min
        if not ma:
            ma = p.max
        if not st:
            st = p.step
        if not se:
            se = p.set
        return np.arange(mi, ma+st, st)

    def path(self, s):
        """Return the path to a file in the results folder"""
        return os.path.join(os.path.dirname(__file__), '../data/results/' + s)

    def w_int(self):

        # Set all values
        inputs = {k: [] for k in self.pars}
        inputs['path'] = []

        # Loop through combinations of min and max pulse widths
        for mi in self.par_range(self.w_int_min):
            for ma in self.par_range(self.w_int_max, mi=mi):

                for p in self.pars:
                    if p == 'w_int_min':
                        inputs[p].append(mi)
                    elif p == 'w_int_max':
                        inputs[p].append(ma)
                    else:
                        inputs[p].append(self.pars[p].set)

                name = 'w_int_{:.2}_{:.2}'.format(mi, ma).replace('.', 'd')
                inputs['path'].append(self.path(name))

        self.dfs['w_int'] = pd.DataFrame(inputs)

    def all_pars(self):
        return pd.concat(self.dfs)

    def save(self, df=None, filename=None):
        """Save database to SQL database"""
        if df is None:
            df = self.all_pars()
        if filename is None:
            filename = 'pars.db'
        conn = sqlite3.connect(self.path(filename))
        df.to_sql('pars', conn, if_exists='replace')
        conn.close()

    def read(self, filename=None):
        """Read in parameter database."""
        if not filename:
            filename = 'pars.db'
        conn = sqlite3.connect(self.path(filename))
        df = pd.read_sql_query("select * from pars;", conn)
        cols = [c for c in df.columns if c.startswith('level')]
        df = df.drop(cols, axis=1)
        return df
