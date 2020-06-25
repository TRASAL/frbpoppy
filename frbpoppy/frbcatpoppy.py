"""Do things with frbcat."""
import os
import numpy as np
import pandas as pd

from frbcat import Frbcat as PureFrbcat

from frbpoppy.misc import pprint
from frbpoppy.paths import paths
from frbpoppy.population import Population


class Frbcat(PureFrbcat):
    """
    Add frbpoppy functionality to Frbcat.

    Get the pandas dataframe with Frbcat().df
    """

    def __init__(self, frbpoppy=True, repeat_bursts=False, **kwargs):
        """Initialize."""
        super().__init__(self, path=paths.frbcat(), **kwargs)

        # Transform the data
        if frbpoppy is True:
            self.frbpoppify()
            self.match_surveys()

        # Just for neating up
        self.df = self.df.sort_values('utc', ascending=False)
        self.df = self.df.reindex(sorted(self.df.columns), axis=1)

    def frbpoppify(self):
        """Prep data for frbpoppy."""
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

        # Gives somewhat of an idea of the pulse width upon arrival at Earth
        self.df['w_arr'] = (self.df['w_eff']**2 -
                            self.df['t_dm']**2 -
                            self.df['t_scat']**2 -
                            self.df['t_samp']**2)**0.5

    def match_surveys(self, interrupt=True):
        """Match up frbs with surveys."""
        # Merge survey names
        surf = os.path.join(self.path, 'paper_survey.csv')
        self._surveys = pd.read_csv(surf)
        cols = ['frb_name', 'pub_description']

        self.df = pd.merge(self.df, self._surveys, on=cols, how='left')
        # Clean up possible unnamed columns
        self.df = self.df.loc[:, ~self.df.columns.str.contains('unnamed')]

        # Add single survey instruments
        # Bit rough, but will work in a pinch
        def cond(t):
            return (self.df.telescope == t) & (self.df.survey.isnull())
        self.df.at[cond('apertif'), 'survey'] = 'apertif'
        self.df.at[cond('askap'), 'survey'] = 'askap-incoh'
        self.df.at[cond('chime'), 'survey'] = 'chime'
        self.df.at[cond('srt'), 'survey'] = 'srt'
        self.df.at[cond('effelsberg'), 'survey'] = 'effelsberg'
        self.df.at[cond('gbt'), 'survey'] = 'guppi'
        self.df.at[cond('fast'), 'survey'] = 'crafts'

        # Check whether any FRBs have not yet been assigned
        no_surveys = self.df['survey'].isnull()

        if interrupt and any(no_surveys):
            cols = ['pub_description', 'frb_name']
            ns_df = self.df[no_surveys].drop_duplicates(subset=cols,
                                                        keep='first')

            pprint('It seems there are new FRBs!')
            m = "  - Frbcat doesn't know which *survey* was running when the "
            m += "FRB was seen"
            pprint(m)
            m = "  - To use these recent detections, please link the FRB to a "
            m += "survey by:"
            pprint(m)
            pprint('  - Adding these frbs to {}'.format(surf))

            for i, r in ns_df[cols].iterrows():
                title, name = r
                if isinstance(title, str):
                    title = title.replace('\n', '')
                print(f'"{title}","{name}",""')

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

        frbs.name = df.name.values
        frbs.dm = df.dm.values
        frbs.dm_mw = df.dm_mw.values
        frbs.gl = df.gl.values
        frbs.gb = df.gb.values
        frbs.ra = df.ra.values
        frbs.dec = df.dec.values
        frbs.z = df.z.values
        frbs.t_scat = df.t_scat.values
        frbs.w_eff = df.w_eff.values
        frbs.si = df.si.values
        frbs.snr = df.snr.values
        frbs.s_peak = df.s_peak.values
        frbs.fluence = df.fluence.values
        # frbs.time = df.utc.values

        pop.save()

        return pop


if __name__ == '__main__':
    import IPython
    f = Frbcat()
    IPython.embed()
