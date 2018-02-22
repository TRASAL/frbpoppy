"""Monte Carlo over various frbpoppy parameters."""

from collections import defaultdict
from scipy.stats import ks_2samp
from datetime import datetime
import numpy as np
import os
import pandas as pd
import sqlite3

from frbpoppy.adapt_pop import Adapt
from frbpoppy.do_hist import histogram
from frbpoppy.do_populate import generate
from frbpoppy.do_survey import observe
from frbpoppy.frbcat import get_frbcat
from frbpoppy.log import pprint


class Parameter:
    """Class to hold attributes of a parameter."""

    def __init__(self, mi, ma, sp, st, log=False):
        """
        Initialisation.

        Args:
            mi (float): Minimum a parameter can be.
            ma (float): Maximum a parameter can be.
            sp (float): Stepsize of the parameter space.
            st (float): Set the default value a parameter should take.
            log (bool, optional): If set to true, it will take the stepsize to
                have been denoted in the power of 10. For instance if the
                minimum has been set to 10**2, and max to 10**5, with step 1,
                then a range in the form [10**2, 10**3, 10**4, 10**5] will be
                produced.
        """
        # Way to check if attributes have changed
        self.altered = -4

        # Set parameters
        self.min = mi
        self.max = ma
        self.step = sp
        self.set = st
        self.log = log

    def __setattr__(self, attrname, val):
        """Check if attributes change after first initialisation."""
        attrs = ['min', 'max', 'step', 'set']
        if attrname in attrs:
            self.altered += 1
        super.__setattr__(self, attrname, val)

    def limits(self, mi, ma, sp, st, log=False):
        """Set parameter range limits."""
        self.min = mi
        self.max = ma
        self.step = sp
        self.set = st
        self.log = log

    def seq(self, start, stop, step):
        """Generate a range of numbers (even if using decimal steps)."""
        n = abs(int(round((stop - start)/float(step))))
        return [start + step*i for i in range(n+1)]

    def par_range(self, mi=None):
        """Quick range generator."""
        # Set lower limit (convenient for double for-loops)
        if not mi:
            mi = self.min

        # Calculate range
        if self.log:
            mi = np.log10(mi)
            ma = np.log10(self.max)
            st = self.step
            r = [10**n for n in self.seq(mi, ma, st)]
        else:
            r = self.seq(mi, self.max, self.step)

        return r


class MonteCarlo:

    def __init__(self):

        # Initialise parameters
        self.dm_host = Parameter(0, 200, 10, 100)
        self.dm_igm_slope = Parameter(1000, 1400, 20, 1200)
        self.freq_max = Parameter(1e6, 1e11, 0.5, 1e10, log=True)
        self.freq_min = Parameter(1e6, 1e11, 0.5, 1e7, log=True)
        self.lum_bol_max = Parameter(1e30, 1e60, 1, 1e50, log=True)
        self.lum_bol_min = Parameter(1e30, 1e60, 1, 1e40, log=True)
        self.lum_bol_slope = Parameter(0.5, 1.5, 0.1, 1.)
        self.n_day = Parameter(2000, 14000, 2000, 10000)
        self.rep = Parameter(0.0, 0.1, 0.01, 0.05)
        self.si_mean = Parameter(-2.0, -1, 0.1, -1.4)
        self.si_sigma = Parameter(0.0, 0.5, 0.1, 0.0)
        self.w_int_max = Parameter(0, 5, 0.1, 5)
        self.w_int_min = Parameter(0, 5, 0.1, 1)

        # Gather parameters
        self.pars = {'dm_host': self.dm_host,
                     'dm_igm_slope': self.dm_igm_slope,
                     'freq_max': self.freq_max,
                     'freq_min': self.freq_min,
                     'lum_bol_min': self.lum_bol_min,
                     'lum_bol_max': self.lum_bol_max,
                     'lum_bol_slope': self.lum_bol_slope,
                     'n_day': self.n_day,
                     'rep': self.rep,
                     'si_mean': self.si_mean,
                     'si_sigma': self.si_sigma,
                     'w_int_max': self.w_int_max,
                     'w_int_min': self.w_int_min,
                     }

        # Set surveys over which to run
        self.surveys = {'WHOLESKY': 'WHOLESKY',
                        'APERTIF': 'APERTIF',
                        'PMSURV': 'parkes',
                        'HTRU': 'parkes',
                        'ASKAP-INCOH': 'ASKAP',
                        'ASKAP-FLY': 'ASKAP',
                        'GBT': 'GBT',
                        'PALFA': 'arecibo',
                        'ARECIBO-SPF': 'arecibo',
                        'ALFABURST': 'arecibo',
                        'UTMOST-1D': 'UTMOST'}

        # Set number of days over which to run the virtual surveys
        self.days = 14

    def path(self, s):
        """Return the path to a file in the results folder."""
        return os.path.join(os.path.dirname(__file__), '../data/results/' + s)

    def save(self, df=None, filename=None, dtypes=None):
        """Save database to SQL database."""
        if df is None:
            df = self.df
        if filename is None:
            filename = 'pars.db'

        conn = sqlite3.connect(self.path(filename))

        # Allow for dtypes to be specified
        if dtypes is None:
            df.to_sql('pars', conn, if_exists='replace')
        else:
            df.to_sql('pars', conn, if_exists='replace', dtype=dtypes)

        # Make index for faster searching
        try:
            index = "CREATE UNIQUE INDEX ix ON pars (id);"
            conn.cursor().execute(index)
        except sqlite3.IntegrityError:
            index = "CREATE INDEX ix ON pars (id);"
            conn.cursor().execute(index)

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

    def possible_pars(self):
        """
        Calculate all possible combinations of parameters.

        Mostly loops over a single parameter while keeping the others at their
        default value, but will do a dual loop over both possible combinations
        of minimum and maximum parameters if both parameters exist.

        Returns:
            df (DataFrame): Table with all parameter combinations.

        """
        inputs = {k: [] for k in self.pars}
        inputs['in_par'] = []

        # For each parameter
        for name, p in self.pars.items():

            # Check whether dual or not
            dual = False
            name_max = name[:-3] + 'max'

            if name.endswith('min'):
                dual = True
            if name.endswith('max'):
                continue

            # Loop over the parameter's range
            for i in p.par_range():

                if dual is False:
                    # Add values
                    inputs[name].append(i)
                    inputs['in_par'].append(name)

                    # Get defaults for all other values
                    for w, v in self.pars.items():
                        if name != w:
                            inputs[w].append(v.set)

                # Set up double loop over limits on a parameter (min and max)
                else:

                    p_max = self.pars[name_max]

                    for j in p_max.par_range(mi=i):

                        # Add values
                        inputs[name].append(i)
                        inputs[name_max].append(j)
                        inputs['in_par'].append(name)

                        # Get defaults for all other values
                        for w, v in self.pars.items():
                            if name != w and name_max != w:
                                inputs[w].append(v.set)

        df = pd.DataFrame(inputs)
        return df

    def run(self):
        """Run a Monte Carlo."""
        # Mercy on anyone trying to understand this code

        # Get the range of parameters over which to loop
        pos_pars = self.possible_pars()
        pos_pars.to_csv('temp.csv')
        # Get the actual observations with which to compare
        cat = get_frbcat()

        # Set up dictionary for results
        d = defaultdict(list)

        # Set up list for binned data
        hists = []

        # Set up dictionary for rates
        rates = defaultdict(list)

        # Iterate over each parameter
        for name, group in pos_pars.groupby('in_par'):

            pprint(name)

            # Generate initial population
            r = group.iloc[0]
            pop = generate(int(r.n_day)*self.days,
                           days=self.days,
                           dm_pars=[r.dm_host, r.dm_igm_slope],
                           emission_pars=[r.freq_min, r.freq_max],
                           lum_dist_pars=[r.lum_bol_min,
                                          r.lum_bol_max,
                                          r.lum_bol_slope],
                           pulse=[r.w_int_min, r.w_int_max],
                           repeat=float(r.rep),
                           si_pars=[r.si_mean, r.si_sigma])

            # Iterate over each value a parameter can take
            for i, r in group.iterrows():

                # Save the values each parameter has been set to
                defaults = r.to_dict()
                # Get an id
                iden = str(datetime.now())
                # Gather survey populations
                sur_pops = []

                # Adapt population
                if name == 'dm_host':
                    pop = Adapt(pop).dm_host(r.dm_host)
                elif name == 'dm_igm_slope':
                    pop = Adapt(pop).dm_igm(r.dm_igm_slope)
                elif name == 'freq_min':
                    pop = Adapt(pop).freq(r.freq_min, r.freq_max)
                elif name == 'lum_bol_min' or name == 'lum_bol_slope':
                    pop = Adapt(pop).lum_bol(r.lum_bol_min,
                                             r.lum_bol_max,
                                             r.lum_bol_slope)
                elif name == 'si_mean' or name == 'si_sigma':
                    pop = Adapt(pop).si(r.si_mean, r.si_sigma)
                elif name == 'w_int_min':
                    pop = Adapt(pop).w_int(r.w_int_min, r.w_int_max)

                for s in self.surveys:
                    # Survey, survey population
                    sur, sur_pop = observe(pop,
                                           s,
                                           output=False,
                                           return_survey=True)

                    if sur_pop.n_srcs == 0:
                        pprint('No {} population'.format(sur_pop.name))
                        continue

                    # Save survey time for rates
                    time = pop.time

                    # Create a pandas dataframe
                    sur_pop = sur_pop.to_df()

                    # Find matching properties
                    cols = [c for c in cat if c in sur_pop]
                    cols.remove('dm_mw')

                    # Collect results
                    ks = {}

                    # KS-test each parameter
                    for c in cols:
                        obs_cat = cat[(cat.survey == s)][c]
                        obs_cat = pd.to_numeric(obs_cat, errors='coerce')
                        obs_pop = pd.to_numeric(sur_pop[c], errors='coerce')

                        ks['ks_' + c] = ks_2samp(obs_pop, obs_cat)[1]

                    # Add as results
                    for p in defaults:
                        d[p].append(defaults[p])
                    for k in ks:
                        d[k].append(ks[k])
                    d['telescope'].append(self.surveys[s])
                    d['survey'].append(s)
                    d['id'].append(iden)

                    # Gather populations to bin
                    if sur_pop.shape[0] > 0:
                        sur_pop['survey'] = s
                        sur_pop['in_par'] = name
                        sur_pop['id'] = iden
                        sur_pop['frbcat'] = False
                        sur_pops.append(sur_pop)

                    # Gather respective frbcat populations to bin
                    cat_pop = cat[(cat.survey == s)].copy()
                    if cat_pop.shape[0] > 0:
                        cat_pop['in_par'] = name
                        cat_pop['id'] = iden
                        cat_pop['frbcat'] = True
                        sur_pops.append(cat_pop)

                    # Gather rates
                    if sur_pop.shape[0] > 0:

                        frbs_per_day = sur.s_frb_rates.det*86400/time
                        srcs_per_day = sur.s_src_rates.det*86400/time

                        rates['survey'].append(s)
                        rates['in_par'].append(name)
                        rates['id'].append(iden)
                        rates['frbs'].append(sur.s_frb_rates.det)
                        rates['sources'].append(sur.s_src_rates.det)
                        rates['days'].append(time/86400)
                        rates['frbs_per_day'].append(frbs_per_day)
                        rates['srcs_per_day'].append(srcs_per_day)

                if sur_pops:
                    hists.append(histogram(sur_pops))

        try:
            db_hists = pd.concat(hists)
        except ValueError:
            raise ValueError('No FRBs seem to have been detected')
        db_ks = pd.DataFrame(d)
        db_rates = pd.DataFrame(rates)

        self.save(df=db_hists, filename='hists_temp5.db')
        self.save(df=db_ks, filename='ks_temp5.db')
        self.save(df=db_rates, filename='rates_temp5.db')
