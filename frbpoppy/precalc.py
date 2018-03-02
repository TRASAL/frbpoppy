"""Create a lookup tables for redshift and the NE2001 dispersion measure."""

import os
import math
import numpy as np
import sqlite3
import sys
from scipy.integrate import quad as integrate

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


def dist_table(dist, H_0=69.6, W_m=0.286, W_v=0.714, z_max=5.0, test=False):
    """
    Create/use a lookup table for distance to redshift.

    Create a list of tuples to lookup the corresponding redshift for a comoving
    distance [Gpc]. Uses formulas from Hoggs et al. (1999) for the cosmological
    calculations, assuming a flat universe. To avoid long calculation times,
    it will check if a previous run with the same parameters has been done,
    which it will then load it. If not, it will calculate a new table, and save
    the table for later runs.

    Args:
        dist (float): Comoving distance [Gpc]
        H_0 (float, optional): Hubble parameter. Defaults to 69.6
        W_m (float, optional): Omega matter. Defaults to 0.286
        W_k (float, optional): Omega vacuum. Defaults to 0.714
        z_max (float, optional): Maximum redshift. Defaults to 5.0
        test (bool): Flag for coarser resolution
    Returns:
        z (float): Redshift

    """
    uni_mods = os.path.join(paths.models(), 'universe/')

    # Initializing
    cl = 299792.458  # Velocity of light [km/sec]

    def cvt(value):
        """Convert a value to a string without a period."""
        return str(value).replace('.', 'd')

    # Filename
    paras = ['h0', cvt(H_0),
             'wm', cvt(W_m),
             'wv', cvt(W_v),
             'zmax', cvt(z_max)]
    f = '-'.join(paras)

    if test:
        step = 0.1
        path = uni_mods + f + '_test.db'
    else:
        step = 0.0001
        path = uni_mods + f + '.db'

    # Setup database
    db = False
    if os.path.exists(path):
        db = True

    # Connect to database
    conn = sqlite3.connect(path)
    c = conn.cursor()

    # Create db
    if not db:

        W_k = 1.0 - W_m - W_v  # Omega curvature

        if W_k != 0.0:
            pprint('Careful - Your cosmological parameters do not sum to 1.0')

        zs = np.arange(0, z_max+step, step)

        # Create database
        c.execute('create table redshift ' +
                  '(dist real, z real)')

        results = []

        # Numerically integrate the following function
        def d_c(x):
            """Comoving distance (Hogg et al, 1999)."""
            return 1/math.sqrt((W_m*(1+x)**3 + W_k*(1+x)**2 + W_v))

        for z in zs:
            d = cl/H_0*integrate(d_c, 0, z)[0]
            d /= 1e3  # Covert from Mpc to Gpc
            results.append((d, z))

            # Give an update on the progress
            sys.stdout.write('\r{}'.format(z))
            sys.stdout.flush()

        # Save results to database
        c.executemany('insert into redshift values (?,?)', results)

        # Make for easier searching
        c.execute('create index ix on redshift (dist)')

        # Save
        conn.commit()

    # Search database
    z = c.execute('select z from redshift where dist > ? limit 1',
                  [dist]).fetchone()[0]

    # Close database
    conn.close()

    return z
