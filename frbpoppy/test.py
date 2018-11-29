import os
import sqlite3
import numpy as np


import frbpoppy.galacticops as go
from frbpoppy import paths


class NE2001Table:
    """Create/use a NE2001 lookup table for dispersion measure."""

    def __init__(self):
        """Initializing."""
        self.set_file_name()

        # Setup database
        self.db = False
        self.rounding = 5
        self.step = 0.1
        self.rounding = 2
        if os.path.exists(self.file_name):
            self.db = True

    def set_file_name(self):
        """Determine filename."""
        uni_mods = os.path.join(paths.models(), 'universe/')
        self.file_name = uni_mods + 'dm_mw.db'

    def create_table(self):
        """Create a lookup table for dispersion meausre"""
        # Connect to database
        conn = sqlite3.connect(self.file_name)
        c = conn.cursor()

        # Set array of coordinates
        gls = np.arange(-180., 180. + self.step, self.step)
        gbs = np.arange(-90., 90. + self.step, self.step)
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

    def lookup(self, gal, gab):
        """Look up associated milky way dispersion measure with gal coords.

        Args:
            gl (array): Galactic longitude [fractional degrees]
            gb (array): Galactic latitude [fractional degrees]

        Returns:
            dm_mw (float): Galactic dispersion measure [pc*cm^-3]

        """
        # Connect to database
        conn = sqlite3.connect(self.file_name)
        c = conn.cursor()

        dm_mw = np.ones_like(gal)

        # Round values
        gal = np.round(gal, self.rounding)
        gab = np.round(gab, self.rounding)

        # Search database
        query = 'select dm_mw from dm where gl=? and gb=? limit 1'

        for i, gl in enumerate(gal):
            dm_mw[i] = c.execute(query, [gl, gab[i]]).fetchone()[0]

        # Close database
        conn.close()

        return dm_mw


if __name__ == '__main__':
    d = CSMDTable()
    d.lookup()
