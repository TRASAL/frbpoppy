"""Do things with frbcat."""

import glob
import io
import os
import pandas as pd
import requests

import frbpoppy.galacticops as go
from frbpoppy.frb import FRB
from frbpoppy.log import pprint
from frbpoppy.paths import paths
from frbpoppy.population import Population
from frbpoppy.source import Source


class Frbcat():
    """
    Do things with the frbcat.

    Mostly used for playing nicely with frbpoppy. Get the pandas dataframe
    with Frbcat().df
    """

    def __init__(self,
                 frbpoppy=True,
                 one_per_frb=False,
                 repeat_bursts=True,
                 repeater=True):
        """Initializing."""
        # Set path
        self.data_dir = paths.frbcat()

        # Get frbcat data
        self.get()

        # Transform the data
        if frbpoppy is True:
            self.filter(one_per_frb=True,
                        repeat_bursts=False,
                        repeater=True)
            self.clean()
            self.coor_trans()
            self.match_surveys()

        # Just for neating up
        self.df = self.df.sort_values('utc', ascending=False)
        self.df = self.df.reindex(sorted(self.df.columns), axis=1)

    def get(self, internet=False, local=True):
        """
        Get frbcat from online or from a local file.

        Args:
            internet (bool): Whether to get a new version of frbcat
            local (bool): Whether to use a local version of frbcat

        """
        if internet:
            # BROKEN AS OF 13/11/2017
            # Frbcat should be updated in January
            # Try using the most recently published frbcat
            try:
                url = 'http://www.astronomy.swin.edu.au/pulsar/frbcat/'
                url += 'table.php?format=text&sep=comma'

                s = requests.get(url).content
                f = io.StringIO(s.decode('utf-8'))

            # Unless there's no internet
            except requests.ConnectionError:
                local = True

        if local:
            # Find latest version of frbcat
            f = max(glob.glob(self.data_dir + '/frbcat*.csv'),
                    key=os.path.getctime)

        self.df = pd.read_csv(f)

    def filter(self,
               one_per_frb=True,
               repeat_bursts=False,
               repeater=True):
        """Filter frbcat in various ways."""
        # Lower all column names
        self.df.columns = map(str.lower, self.df.columns)

        if one_per_frb:
            # Only keep rows with the largest number of parameters
            # so that only one row per FRB remains
            self.df['count'] = self.df.count(axis=1)
            self.df = self.df.sort_values('count', ascending=False)
            self.df = self.df.drop_duplicates(subset=['utc'])

        if not repeat_bursts:
            # Only keeps one detection of the repeater
            self.df = self.df.drop_duplicates(subset=['frb'])

        if not repeater:
            self.df = self.df[self.df.frb_name != 'FRB121102']

        self.df = self.df.sort_index()

    def clean(self):
        """Clean up the data."""
        # Clean up column names
        self.df.columns = self.df.columns.str.replace('rop_', '')
        self.df.columns = self.df.columns.str.replace('rmp_', '')

        # Conversion table
        convert = {'mw_dm_limit': 'dm_mw',
                   'width': 'w_eff',
                   'flux': 's_peak',
                   'redshift_host': 'z',
                   'spectral_index': 'si',
                   'dispersion_smearing': 't_dm',
                   'dm_error': 'dm_err',
                   'scattering_timescale': 't_scat',
                   'sampling_time': 't_samp'}

        self.df.rename(columns=convert, inplace=True)

        # Add some extra columns
        self.df['fluence'] = self.df['s_peak'] * self.df['w_eff']
        self.df['population'] = 'frbcat'
        self.df['t_dm_err'] = ((self.df['t_dm']/self.df['dm']) *
                               (self.df['dm_err']*self.df['dm']))

        # Temporary fix to missing data on FRB170827 in frbcat
        self.df['t_samp'] = self.df['t_samp'].fillna(0.06400)

        # Gives somewhat of an idea of the pulse width upon arrival at Earth
        self.df['w_arr'] = (self.df['w_eff']**2 -
                            self.df['t_dm']**2 -
                            self.df['t_dm_err']**2 -
                            self.df['t_scat']**2 -
                            self.df['t_samp']**2)**0.5

        # Reduce confusion in telescope names
        small_tele = self.df['telescope'].str.lower()
        self.df['telescope'] = small_tele

        # Split out errors on values
        for c in self.df.columns:
            if self.df[c].dtype == object:
                if any(self.df[c].str.contains('&plusmn', na=False)):
                    val, err = self.df[c].str.split('&plusmn', 1).str
                    self.df[c] = pd.to_numeric(val)
                    self.df[c+'_err'] = pd.to_numeric(err)

        # Set utc as dates
        self.df['utc'] = pd.to_datetime(self.df['utc'])

    def coor_trans(self):
        """Apply coordinate transformations."""
        def trans(df):

            # Clean up some errors in frbcat
            if df['decj'].count(':') < 2:
                df['decj'] = df['decj'] + ':00'

            ra, dec = go.frac_deg(df['raj'], df['decj'])
            gl, gb = go.radec_to_lb(ra, dec, frac=True)
            df['ra'] = ra
            df['dec'] = dec
            df['gl'] = gl
            df['gb'] = gb
            return df

        self.df = self.df.apply(trans, axis=1)

    def match_surveys(self):
        """Match up frbs with surveys."""
        # Merge survey names
        surf = os.path.join(self.data_dir, 'frb_survey.csv')
        self._surveys = pd.read_csv(surf)
        self.df = pd.merge(self.df, self._surveys, on='frb', how='outer')

        # Check whether any FRBs have not yet been assigned
        no_surveys = self.df['survey'].isnull()

        if any(no_surveys):
            names = []
            for index, row in self.df[no_surveys].iterrows():
                names.append(row['frb'])
                m = '====> {} has not been linked to a survey <===='
                n = m.format(row['frb'])
                pprint(n)
                pprint(row)

            pprint('Please add these frbs to {}'.format(surf))

            for name in names:
                pprint(name + ',')

    def to_pop(self, df=None):
        """
        Convert to a Population object.

        Please ensure self.clean() has been run first.
        """
        if not isinstance(df, pd.DataFrame):
            df = self.df

        pop = Population()
        pop.name = 'frbcat'

        # For each source
        for name, src_df in df.groupby('frb'):

            source = Source()
            source.name = name
            source.dm = src_df.dm.iloc[0]
            source.dm_mw = src_df.dm_mw.iloc[0]
            source.gl = src_df.gl.iloc[0]
            source.gb = src_df.gb.iloc[0]
            source.ra = src_df.ra.iloc[0]
            source.dec = src_df.dec.iloc[0]
            source.z = src_df.z.iloc[0]
            source.t_scat = src_df.t_scat.iloc[0]

            for index, row in src_df.iterrows():
                frb = FRB()
                frb.w_eff = row.w_eff
                frb.si = row.si
                frb.snr = row.snr
                frb.s_peak = row.s_peak
                frb.fluence = row.fluence
                source.add(frb)

            pop.add(source)

        pop.save()

        return pop
