"""Do things with frbcat."""

import glob
import io
import os
import pandas as pd
import requests

import frbpoppy.galacticops as go
from frbpoppy.log import pprint
from frbpoppy.paths import paths


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
        if frbpoppy:
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
        if one_per_frb:
            # Only keep rows with the largest number of parameters
            # so that only one row per FRB remains
            self.df['count'] = self.df.count(axis=1)
            self.df = self.df.sort_values('count', ascending=False)
            self.df = self.df.drop_duplicates('utc')

        if not repeat_bursts:
            # Only keeps one detection of the repeater
            self.df = self.df.drop_duplicates('frb')

        if not repeater:
            self.df = self.df[self.df.frb_name != 'FRB121102']

        self.df = self.df.sort_index()

    def clean(self):
        """Clean up the data."""
        # Clean up column names
        self.df.columns = self.df.columns.str.replace('rop_', '')
        self.df.columns = self.df.columns.str.replace('rmp_', '')

        # Conversion table
        convert = {'frb_name': 'frb',
                   'mw_dm_limit': 'dm_mw',
                   'width': 'w_eff',
                   'flux': 's_peak',
                   'redshift_host': 'z',
                   'spectral_index': 'si',
                   'dispersion smearing': 't_dm',
                   'scattering_timescale': 't_scat'}

        self.df.rename(columns=convert, inplace=True)

        # Add some extra columns
        self.df['fluence'] = self.df['s_peak'] * self.df['w_eff']
        self.df['population'] = 'frbcat'

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
