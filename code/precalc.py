"""Create a lookup table for various galactic operations."""

import os
import numpy as np
import sqlite3
import sys

import galacticops as go


mods = os.path.join(os.path.dirname(__file__), '../data/models/')
uni_mods = os.path.join(mods, 'universe/')
dm_mods = os.path.join(mods, 'dm/')


def ne2001_table(gal, gab):
    """
    Create/use a NE2001 lookup table for dispersion measure.

    Args:
        gl (float): Galactic longitude [fractional degrees]
        gb (float): Galactic latitude [fractional degrees]
        df (DataFrame): Lookup table (see coords_table())

    Returns:
        dm_mw (float): Galactic dispersion measure [pc*cm^-3]

    """
    # Setup database
    db = False
    path = uni_mods + 'dm_mw.db'
    if os.path.exists(path):
        db = True

    # Connect to database
    conn = sqlite3.connect(path)
    c = conn.cursor()

    # Create db
    if not db:
        # Set array of coordinates
        step = 0.1
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
    def frac_round(x, prec=2, base=1):
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


ne2001_table()
