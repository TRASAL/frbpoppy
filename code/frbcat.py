"""Import frbcat."""

from numpy import loadtxt
import io
import glob
import os
import pandas as pd
import requests

import galacticops as go


def transform_coords(db):
    """Apply coordinate transformations to frbcat."""
    # Fixing error in frbcat
    if db['DECJ'].count(':') < 2:
        db['DECJ'] = db['DECJ'].replace('.', ':')
        if db['DECJ'].count(':') < 2:
            db['DECJ'] = db['DECJ'] + ':00'

    gl, gb = go.radec_to_lb(db['RAJ'], db['DECJ'])

    db['gl'] = gl
    db['gb'] = gb

    ra, dec = go.lb_to_radec(db['gl'], db['gb'])

    db['ra'] = ra
    db['dec'] = dec

    return db


def match_surveys(db, surveys):
    """Match frbcat entries with a survey."""
    db['survey'] = surveys[db['Name']]
    return db


def get_frbcat():
    """Get frbcat as a Pandas DataFrame."""
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
        f = min(glob.glob(catdir + '/frbcat*.csv'), key=os.path.getctime)

    db = pd.read_csv(f)

    # Only keep rows with the largest number of parameters
    # so that only one row per FRB remains
    db['count'] = db.count(axis=1)
    db = db.sort_values('count', ascending=False).drop_duplicates('UTC')
    db = db.sort_index()

    # Apply coordinate transformations
    db = db.apply(transform_coords, axis=1)

    # Change some names
    db['dm'] = db['DM']
    db['dm_mw'] = db['NE2001 DM Limit']
    db['snr'] = db['SNR']
    db['w_eff'] = db['Width']
    db['s_peak'] = db['Flux']
    db['fluence'] = db['Flux'] * db['Width']
    db['population'] = 'frbcat'
    db['z'] = db['Spectroscopic Redshift']

    # Match up with surveys
    cwd = os.path.dirname(__file__)
    folder = '../data/frbcat/'
    surdir = os.path.join(cwd, folder, 'frb_survey.csv')
    key_value = loadtxt(surdir, dtype=str, delimiter=",")
    surveys = {k: v for k, v in key_value}
    db = db.apply(match_surveys, args=(surveys,), axis=1)

    return db
