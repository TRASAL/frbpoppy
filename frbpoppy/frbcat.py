"""Do things with frbcat."""
import datetime
import glob
import io
import os
import pandas as pd
import requests

from frbpoppy.log import pprint
from frbpoppy.paths import paths
from frbpoppy.population import Population
import frbpoppy.galacticops as go


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
                 repeaters=True,
                 update=False):
        """Initializing."""
        # Set path
        self.data_dir = paths.frbcat()

        # Get frbcat data
        self.get(internet=update)

        # Transform the data
        if frbpoppy is True:
            self.clean()
            self.filter(one_per_frb=True,
                        repeat_bursts=False,
                        repeaters=True)
            self.coor_trans()
            self.match_surveys()

        # Just for neating up
        self.df = self.df.sort_values('utc', ascending=False)
        self.df = self.df.reindex(sorted(self.df.columns), axis=1)

    def url_to_df(self, url):
        """Convert a url of a JSON table to a Pandas DataFrame.

        Args:
            url (str): URL to the webpage

        Returns:
            DataFrame: DataFrame of JSON table

        """
        try:
            s = requests.get(url).content
            f = io.StringIO(s.decode('utf-8'))

            series = []
            for entry in pd.read_json(f)['products']:
                series.append(pd.Series(entry))
            df = pd.concat(series, axis=1).T

            return df

        except ValueError:
            pass

    def urls_to_df(self, endings, url):
        """
        Use Series to loop over multiple webpages and concatenate them to a
        single DataFrame

        Args:
            endings (iterables): The list/series/column over which to loop
            url (str): The base url

        Returns:
            DataFrame

        """
        dfs = []
        for ending in endings:
            full_url = f'{url}{ending}'
            df = self.url_to_df(full_url)
            dfs.append(df)
        return pd.concat(dfs, ignore_index=True)

    def get(self, internet=False, save=True, local=False):
        """
        Get frbcat from online or from a local file.

        Args:
            internet (bool): Whether to get a new version of frbcat
            local (bool): Whether to use a local version of frbcat

        """
        # Check whether a copy of FRBCAT has already been downloaded
        # Ensures frbcat is only queried once a month
        path = self.data_dir + '/frbcat_'
        path += str(datetime.datetime.today()).split()[0][:-3]
        path += '-??.csv'
        exists = glob.glob(path)
        if internet and exists:
            internet = False
            local = True

        if internet:
            try:
                pprint('Attempting to retrieve FRBCAT from www.frbcat.org')

                # First get all FRB names from the main page
                url = 'http://frbcat.org/products/'
                main_df = self.url_to_df(url)

                # Then get any subsequent analyses (multiple entries per FRB)
                url = 'http://frbcat.org/product/'
                frb_df = self.urls_to_df(main_df.frb_name, url)

                # Find all frb note properties
                url = 'http://frbcat.org/frbnotes/'
                frbnotes_df = self.urls_to_df(set(frb_df.index), url)
                frbnotes_df = frbnotes_df.add_prefix('frb_notes_')

                # Find all notes on radio observation parameters
                url = 'http://frbcat.org/ropnotes/'
                ropnotes_df = self.urls_to_df(set(frb_df.index), url)
                ropnotes_df = ropnotes_df.add_prefix('rop_notes_')

                # Find all radio measurement parameters
                url = 'http://frbcat.org/rmppubs/'
                rmppubs_df = self.urls_to_df(set(frb_df.index), url)
                rmppubs_df = rmppubs_df.add_prefix('rmp_pub_')

                # Have skipped
                # 'http://frbcat.org/rmpimages/<rmp_id>' (images)
                # 'http://frbcat.org/rmpnotes/<rmp_id>' (empty)

                # Merge all databases together
                df = pd.merge(frb_df,
                              frbnotes_df,
                              left_on='frb_id',
                              right_on='frb_notes_frb_id',
                              how='left')

                df = pd.merge(df,
                              ropnotes_df,
                              left_on='rop_id',
                              right_on='rop_notes_rop_id',
                              how='left')

                self.df = pd.merge(df,
                                   rmppubs_df,
                                   left_on='rmp_id',
                                   right_on='rmp_pub_rmp_id',
                                   how='left')

                pprint('Succeeded')

                if save:
                    date = str(datetime.datetime.today()).split()[0]
                    path = self.data_dir + f'/frbcat_{date}.csv'
                    self.df.to_csv(path)

                local = False

            # Unless there's no internet
            except requests.ConnectionError:
                local = True

        if local or not internet:
            # Find latest version of frbcat
            f = max(glob.glob(self.data_dir + '/frbcat*.csv'),
                    key=os.path.getctime)
            print(f)
            self.df = pd.read_csv(f)

    def clean(self):
        """Clean up the data."""
        # Lower all column names
        self.df.columns = map(str.lower, self.df.columns)

        # Convert None's to Nan's
        self.df.fillna(value=pd.np.nan, inplace=True)

        # Clean up column names
        self.df.columns = self.df.columns.str.replace('rop_', '')
        self.df.columns = self.df.columns.str.replace('rmp_', '')

        # There's a problem with mulitple 'id' columns
        cols = [c for c in self.df.columns if not c.endswith('id')]
        self.df = self.df[cols]

        # Split out errors on values
        for c in self.df.columns:
            if self.df[c].dtype == object:
                if any(self.df[c].str.contains('&plusmn', na=False)):
                    val, err = self.df[c].str.split('&plusmn', 1).str
                    self.df[c] = pd.to_numeric(val)
                    self.df[c+'_err'] = pd.to_numeric(err)

        # Split out asymetric errors on values
        for c in self.df.columns:
            if self.df[c].dtype == object:
                if any(self.df[c].str.contains('<sup>', na=False)):
                    upper = "<span className='supsub'><sup>"
                    val, rest = self.df[c].str.split(upper, 1).str
                    upper, rest = rest.str.split('</sup><sub>', 1).str
                    lower, _ = rest.str.split('</sub></span>', 1).str
                    self.df[c] = pd.to_numeric(val)
                    self.df[c+'_err_up'] = pd.to_numeric(upper)
                    self.df[c+'_err_down'] = pd.to_numeric(lower)

        # Conversion table
        convert = {'frb_name': 'frb',
                   'mw_dm_limit': 'dm_mw',
                   'width': 'w_eff',
                   'flux': 's_peak',
                   'redshift_host': 'z',
                   'spectral_index': 'si',
                   'dispersion_smearing': 't_dm',
                   'dm_error': 'dm_err',
                   'scattering_timescale': 't_scat',
                   'sampling_time': 't_samp'}

        self.df.rename(columns=convert, inplace=True)

        # Ensure columns are the right datatype
        self.df.w_eff = pd.to_numeric(self.df.w_eff, errors='coerce')

        # Add some extra columns
        self.df['fluence'] = self.df['s_peak'] * self.df['w_eff']
        self.df['population'] = 'frbcat'

        # Gives somewhat of an idea of the pulse width upon arrival at Earth
        self.df['w_arr'] = (self.df['w_eff']**2 -
                            self.df['t_dm']**2 -
                            self.df['t_scat']**2 -
                            self.df['t_samp']**2)**0.5

        # Reduce confusion in telescope names
        small_tele = self.df['telescope'].str.lower()
        self.df['telescope'] = small_tele

        # Set utc as dates
        self.df['utc'] = pd.to_datetime(self.df['utc'])

        # Replace chime/frb with chime
        if any(self.df['telescope'].str.contains('chime/frb', na=False)):
            val, _ = self.df['telescope'].str.split('/', 1).str
            self.df['telescope'] = val

    def filter(self,
               one_per_frb=True,
               repeat_bursts=False,
               repeaters=True):
        """Filter frbcat in various ways."""
        if one_per_frb:
            # Only keep rows with the largest number of parameters
            # so that only one row per detected FRB remains
            self.df['count'] = self.df.count(axis=1)
            self.df = self.df.sort_values('count', ascending=False)
            self.df = self.df.drop_duplicates(subset=['utc'])

        if not repeaters:
            # Drops any repeater sources
            self.df = self.df.drop_duplicates(subset=['frb'], keep=False)

        if not repeat_bursts:
            # Only keeps one detection of repeaters
            self.df = self.df.drop_duplicates(subset=['frb'])

        self.df = self.df.sort_index()

    def coor_trans(self):
        """Apply coordinate transformations."""
        def trans(df):

            # Clean up some errors in frbcat
            if df['decj'].count(':') < 2:
                df['decj'] = df['decj'] + ':00'
            if df['raj'].count(':') < 2:
                df['raj'] = df['raj'] + ':00'

            ra, dec = go.frac_deg(df['raj'], df['decj'])
            gl, gb = go.radec_to_lb(ra, dec, frac=True)
            df['ra'] = ra
            df['dec'] = dec
            df['gl'] = gl
            df['gb'] = gb
            return df

        self.df = self.df.apply(trans, axis=1)

    def match_surveys(self, interrupt=True):
        """Match up frbs with surveys."""
        # Merge survey names
        surf = os.path.join(self.data_dir, 'paper_survey.csv')
        self._surveys = pd.read_csv(surf)
        cols = ['frb', 'pub_description']
        self.df = pd.merge(self.df, self._surveys, on=cols, how='outer')

        # Clean up possible unnamed columns
        self.df = self.df.loc[:, ~self.df.columns.str.contains('unnamed')]

        # Check whether any FRBs have not yet been assigned
        no_surveys = self.df['survey'].isnull()

        if interrupt:
            if any(no_surveys):
                ns_df = self.df[no_surveys]

                pprint('Please add these frbs to {}'.format(surf))

                for i, r in ns_df[['pub_description', 'frb']].iterrows():
                    title, frb = r
                    print(f'"{title}",{frb},')

    def to_pop(self, df=None):
        """
        Convert to a Population object.

        Please ensure self.clean() has been run first.
        """
        if not isinstance(df, pd.DataFrame):
            df = self.df

        pop = Population()
        pop.name = 'frbcat'
        frbs = pop.frbs

        frbs.name = df.frb
        frbs.dm = df.dm
        frbs.dm_mw = df.dm_mw
        frbs.gl = df.gl
        frbs.gb = df.gb
        frbs.ra = df.ra
        frbs.dec = df.dec
        frbs.z = df.z
        frbs.t_scat = df.t_scat
        frbs.w_eff = df.w_eff
        frbs.si = df.si
        frbs.snr = df.snr
        frbs.s_peak = df.s_peak
        frbs.fluence = df.fluence

        pop.save()

        return pop


if __name__ == '__main__':
    f = Frbcat()
