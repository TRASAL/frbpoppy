import io
import glob
import os
import pandas as pd
import requests

import galacticops as go


def transform_coords(db):
    """Apply coordinate transformations to frbcat"""

    # Fixing error in frbcat
    if db['DECJ'].count(':') < 2:
        db['DECJ'] = db['DECJ'].replace('.', ':')
        if db['DECJ'].count(':') < 2:
            db['DECJ'] = db['DECJ'] + ':00'

    gl, gb = go.radec_to_lb(db['RAJ'], db['DECJ'])

    db['gl'] = gl
    db['gb'] = gb

    return db


def get_frbcat():

    # Try using the most recently published frbcat
    try:
        url = 'http://www.astronomy.swin.edu.au/pulsar/frbcat/'
        url += 'table.php?format=text&sep=comma'

        s = requests.get(url).content
        f = io.StringIO(s.decode('utf-8'))

    # Unless there's no internet
    except requests.ConnectionError:
        cwd = os.path.dirname(__file__)
        folder = '../data/frbcat/'
        catdir = os.path.join(cwd, folder)
        # Find latest version of frbcat
        f = min(glob.glob(catdir + '/*.csv'), key=os.path.getctime)

    db = pd.read_csv(f)

    # Only keep rows with the largest number of parameters
    db['count'] = db.count(axis=1)
    db = db.sort_values('count', ascending=False).drop_duplicates('UTC')
    db = db.sort_index()

    # Apply coordinate transformations
    db = db.apply(transform_coords, axis=1)

    # Change some names
    db['dm'] = db['DM']
    db['dm_mw'] = db['NE2001 DM Limit']
    db['si'] = db['Scattering Index']
    db['snr'] = db['SNR']
    db['w_eff'] = db['Width']
    db['s_peak'] = 1e26*db['Flux']  # Convert to W/(m^2*Hz)
    db['fluence'] = db['Flux'] * db['Width']
    db['population'] = 'frbcat'
    db['z'] = db['Spectroscopic Redshift']

    return db


if __name__ == '__main__':
    get_frbcat()
