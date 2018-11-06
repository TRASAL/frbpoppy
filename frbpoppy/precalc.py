"""Create a lookup tables for redshift and the NE2001 dispersion measure."""

import os
import numpy as np
import sqlite3
import sys
from scipy.integrate import quad

import frbpoppy.galacticops as go
from frbpoppy.log import pprint
from frbpoppy.paths import paths


def ne2001_table(gal, gab, test=False):
    """
    Create/use a NE2001 lookup table for dispersion measure.

    Args:
        gl (float): Galactic longitude [fractional degrees]
        gb (float): Galactic latitude [fractional degrees]
        test (bool): Flag for coarser resolution

    Returns:
        dm_mw (float): Galactic dispersion measure [pc*cm^-3]

    """
    uni_mods = os.path.join(paths.models(), 'universe/')

    # Set up for testing
    if test:
        step = 10
        rounding = -1
        path = uni_mods + 'dm_mw_test.db'
    else:
        step = 0.1
        rounding = 2
        path = uni_mods + 'dm_mw.db'

    # Setup database
    db = False
    if os.path.exists(path):
        db = True

    # Connect to database
    conn = sqlite3.connect(path)
    c = conn.cursor()

    # Create db
    if not db:
        # Set array of coordinates
        gls = np.arange(-180., 180. + step, step)
        gbs = np.arange(-90., 90. + step, step)
        dist = 0.1  # [Gpc]

        # Create database
        c.execute('create table dm ' +
                  '(gl real, gb real, dm_mw real)')

        results = []

        for gl in gls:
            gl = round(gl, 1)
            for gb in gbs:
                gb = round(gb, 1)

                dm_mw = go.ne2001_dist_to_dm(dist, gl, gb)

                r = (gl, gb, dm_mw)
                results.append(r)

            # Give an update on the progress
            pprint('Creating a DM lookup table (only needs to be done once)')
            sys.stdout.write('\r{}'.format(gl))
            sys.stdout.flush()

        # Save results to database
        c.executemany('insert into dm values (?,?,?)', results)

        # Make for easier searching
        c.execute('create index ix on dm (gl, gb)')

        # Save
        conn.commit()

    # Round values
    def frac_round(x, prec=rounding, base=1):
        return round(base * round(float(x)/base), prec)

    # Round to 0.05 fractional degree
    gal = frac_round(gal)
    gab = frac_round(gab)

    # Search database
    dm_mw = c.execute('select dm_mw from dm where gl=? and gb=? limit 1',
                      [gal, gab]).fetchone()[0]

    # Close database
    conn.close()

    return dm_mw


def dist_table(H_0=69.6, W_m=0.286, W_v=0.714, test=False,
               dist_co=None, vol_co=None, z=None):
    """
    Create/use a lookup table for comoving distance, volume & redshift.

    Create a list of tuples to lookup the corresponding redshift for a comoving
    distance [Gpc] (or the other way around). Uses formulas from
    Hoggs et al. (1999) for the cosmological calculations. To avoid long
    calculation times, it will check if a previous run with the same parameters
    has been done, which it will then load it. If not, it will calculate a new
    table, and save the table for later runs.

    Args:
        H_0 (float, optional): Hubble parameter. Defaults to 69.6 km/s/Mpc
        W_m (float, optional): Omega matter. Defaults to 0.286
        W_k (float, optional): Omega vacuum. Defaults to 0.714
        test (bool): Flag for coarser resolution
        dist_co (float or bool): Whether to give or return comoving distance
        z (float or bool): Whether to give or return redshift
        vol_co (float or bool): Whether to give or comoving volume

    Returns:
        float: Distance measures [Gpc], redshift, or comoving volume from
            Earth to the source [Gpc^3]

        Alternatively
        dict: Various outputs as requested

    """
    uni_mods = os.path.join(paths.models(), 'universe/')

    # Check input
    if not any([e for e in [dist_co, vol_co, z] if e is True]):
        raise ValueError('Set which parameters you would like as output')
    # Using z later on
    redshift = z

    def cvt(value):
        """Convert a float to a string without a period."""
        return str(value).replace('.', 'd')

    # Filename
    paras = ['h0', cvt(H_0),
             'wm', cvt(W_m),
             'wv', cvt(W_v)]
    f = '-'.join(paras)

    if test:
        step = 0.1
        z_max = 3.0
        file_name = uni_mods + f + '_test.db'
    else:
        step = 0.0001
        z_max = 8.0
        file_name = uni_mods + f + '.db'

    # Setup database
    db = False
    if os.path.exists(file_name):
        db = True

    # Connect to database
    conn = sqlite3.connect(file_name)
    c = conn.cursor()

    # Create db
    if not db:

        pprint('Creating a distance table (only needs to happen once)')
        pprint('At redshift:')

        W_k = 1.0 - W_m - W_v  # Omega curvature

        if W_k != 0.0:
            pprint('Careful - Your cosmological parameters do not sum to 1.0')

        zs = np.arange(0, z_max+step, step)

        # Create database
        c.execute('create table distances (z real, dist real, vol real)')

        results = []

        for z in zs:
            # Comoving distance [Gpc]
            dist = go.Redshift(z, H_0=H_0, W_m=W_m, W_v=W_v)
            results.append((z, dist.dist_co(), dist.vol_co()))

            # Give an update on the progress
            sys.stdout.write('\r{}'.format(z))
            sys.stdout.flush()

        # Save results to database
        c.executemany('insert into distances values (?,?,?)', results)

        # Make for easier searching
        c.execute('create index ix on distances (z)')
        c.execute('create index ixx on distances (dist)')
        c.execute('create index ixxx on distances (vol)')

        # Save
        conn.commit()

        pprint('\nFinished distance table')

    # Set input value
    if isinstance(dist_co, float):
        in_par = 'dist'
        i = dist_co
    elif isinstance(vol_co, float):
        in_par = 'vol'
        i = vol_co
    elif isinstance(redshift, float):
        in_par = 'z'
        i = redshift

    # Search database
    query = f'select * from distances where {in_par} > ? limit 1'
    results = c.execute(query, [i]).fetchone()

    # Close database
    conn.close()

    # Gather outputs
    outputs = {}
    if redshift is True:
        outputs['z'] = results[0]
    if dist_co is True:
        outputs['dist_co'] = results[1]
    if vol_co is True:
        outputs['vol_co'] = results[2]

    if len(outputs) == 1:
        output, = outputs.values()
        return output
    else:
        return outputs


def csmd_table(z, H_0=69.6, W_m=0.286, W_v=0.714, test=False):
    """
    Create/use a stellar mass density lookup table for the FRB number density.

    Args:
        z (float): Redshift
        H_0 (float, optional): Hubble parameter. Defaults to 69.6
        W_m (float, optional): Omega matter. Defaults to 0.286
        W_v (float, optional): Omega vacuum. Defaults to 0.714
        test (bool): Flag for coarser resolution

    Returns:
        float: Stellar mass density at given redshift

    """
    uni_mods = os.path.join(paths.models(), 'universe/')

    # Set up for testing
    if test:
        step = 0.1
        rounding = 1
        path = uni_mods + 'csmd_test.db'
    else:
        step = 1e-5
        rounding = 5
        path = uni_mods + 'csmd.db'

    # Setup database
    db = False
    if os.path.exists(path):
        db = True

    # Connect to database
    conn = sqlite3.connect(path)
    c = conn.cursor()

    # Create db
    if not db:
        # Set array of coordinates
        zs = np.arange(0., 6. + step, step)

        # Create database
        c.execute('create table csmd (z real, csmd real)')

        results = []

        def integral(z):
            z1 = z + 1
            return z1**1.7/(1+(z1/2.9)**5.6)*(1/(H_0*(W_m*z1**3+W_v)**0.5))

        def csmd(z):
            return 0.01095*quad(integral, z, np.inf)[0]

        vec_csmd = np.vectorize(csmd)

        # Give an update on the progress
        pprint('Creating a CSMD lookup table (only needs to be done once)')

        csmds = vec_csmd(zs)

        results = list(zip(np.around(zs, decimals=rounding), csmds))

        # Save results to database
        c.executemany('insert into csmd values (?,?)', results)

        # Make for easier searching
        c.execute('create index ix on csmd (z, csmd)')

        # Save
        conn.commit()

    # Round values
    z = round(z, rounding)

    # Search database
    smd = c.execute('select csmd from csmd where z=? limit 1',
                    [z]).fetchone()[0]

    # Close database
    conn.close()

    return smd
