"""Create a lookup tables for redshift and the NE2001 dispersion measure."""

import os
import numpy as np
import sqlite3
import sys
from scipy.integrate import quad

import frbpoppy.galacticops as go
from frbpoppy.log import pprint
from frbpoppy.paths import paths


class NE2001Table:
    """Create/use a NE2001 lookup table for dispersion measure."""

    def __init__(self):
        """Initializing."""
        self.set_file_name()

        # Setup database
        self.db = False
        self.step = 0.1
        self.rounding = 2
        if os.path.exists(self.file_name):
            self.db = True
        else:
            self.create_table()

    def set_file_name(self):
        """Determine filename."""
        uni_mods = os.path.join(paths.models(), 'universe/')
        self.file_name = uni_mods + 'dm_mw.db'

    def create_table(self):
        """Create a lookup table for dispersion measure."""
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

        # Give an update on the progress
        pprint('Creating a DM lookup table (only needs to be done once)')

        for gl in gls:
            gl = round(gl, 1)
            for gb in gbs:
                gb = round(gb, 1)

                dm_mw = go.ne2001_dist_to_dm(dist, gl, gb)

                r = (gl, gb, dm_mw)
                results.append(r)

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
        def frac_round(x, prec=self.rounding, base=1):
            return np.round(base * np.round(x/base), prec)

        # Round values
        gal = frac_round(gal, self.rounding)
        gab = frac_round(gab, self.rounding)

        # Search database
        query = 'select dm_mw from dm where gl=? and gb=? limit 1'

        for i, gl in enumerate(gal):
            dm_mw[i] = c.execute(query, [gl, gab[i]]).fetchone()[0]

        # Close database
        conn.close()

        return dm_mw


class DistanceTable:
    """
    Create/use a lookup table for comoving distance, volume & redshift.

    Create a list of tuples to lookup the corresponding redshift for a comoving
    distance [Gpc] (or the other way around). Uses formulas from
    Hoggs et al. (1999) for the cosmological calculations. To avoid long
    calculation times, it will check if a previous run with the same parameters
    has been done, which it will then load it. If not, it will calculate a new
    table, and save the table for later runs.

    Args:
        H_0 (float, optional): Hubble parameter. Defaults to 67.74 km/s/Mpc
        W_m (float, optional): Omega matter. Defaults to 0.3089
        W_k (float, optional): Omega vacuum. Defaults to 0.6911

    """

    def __init__(self, H_0=67.74, W_m=0.3089, W_v=0.6911):
        """Initializing."""
        self.H_0 = H_0
        self.W_m = W_m
        self.W_v = W_v

        self.set_file_name()

        # Setup database
        self.db = False
        self.step = 0.00001
        self.z_max = 6.5
        if os.path.exists(self.file_name):
            self.db = True
        else:
            self.create_table()

    def set_file_name(self):
        """Determine filename."""
        uni_mods = os.path.join(paths.models(), 'universe/')

        def cvt(value):
            """Convert a float to a string without a period."""
            return str(value).replace('.', 'd')

        # Convert
        paras = ['h0', cvt(self.H_0),
                 'wm', cvt(self.W_m),
                 'wv', cvt(self.W_v)]
        f = '-'.join(paras)

        self.file_name = uni_mods + f + '.db'

    def create_table(self):
        """Create a lookup table for distances."""
        pprint('Creating a distance table (only needs to happen once)')

        # Connect to database
        conn = sqlite3.connect(self.file_name)
        c = conn.cursor()

        H_0 = self.H_0
        W_m = self.W_m
        W_v = self.W_v

        W_k = 1.0 - W_m - W_v  # Omega curvature

        if W_k != 0.0:
            pprint('Careful - Your cosmological parameters do not sum to 1.0')

        zs = np.arange(0, self.z_max+self.step, self.step)

        # Create database
        c.execute('create table distances (z real, dist real, vol real)')

        results = []

        pprint('At redshift:')

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

    def lookup(self, z=None, dist_co=None, vol_co=None):
        """Look up associated values with input values."""
        # Connect to database
        conn = sqlite3.connect(self.file_name)
        c = conn.cursor()

        # Check what's being looked up
        if z is not None:
            in_par = 'z'
            dist_co = np.ones_like(z)
            vol_co = np.ones_like(z)
        elif dist_co is not None:
            in_par = 'dist'
            z = np.ones_like(dist_co)
            vol_co = np.ones_like(dist_co)
        else:
            in_par = 'vol'
            z = np.ones_like(vol_co)
            dist_co = np.ones_like(vol_co)

        # Search database
        query = f'select * from distances where {in_par} > ? limit 1'

        if in_par == 'z':
            for i, r in enumerate(z):
                _, dist_co[i], vol_co[i] = c.execute(query, [r]).fetchone()
        if in_par == 'dist':
            for i, r in enumerate(dist_co):
                z[i], _, vol_co[i] = c.execute(query, [r]).fetchone()
        if in_par == 'vol':
            for i, r in enumerate(vol_co):
                z[i], dist_co[i], _ = c.execute(query, [r]).fetchone()

        # Close database
        conn.close()

        return z, dist_co, vol_co


class CSMDTable:
    """
    Create/use a stellar mass density lookup table for the FRB number density.

    Args:
        H_0 (float, optional): Hubble parameter. Defaults to 67.74
        W_m (float, optional): Omega matter. Defaults to 0.3089
        W_v (float, optional): Omega vacuum. Defaults to 0.6911

    """

    def __init__(self, H_0=67.74, W_m=0.3089, W_v=0.6911):
        """Initializing."""
        self.H_0 = H_0
        self.W_m = W_m
        self.W_v = W_v

        self.set_file_name()

        # Setup database
        self.db = False
        self.rounding = 5
        self.step = 1e-5
        if os.path.exists(self.file_name):
            self.db = True
        else:
            self.create_table()

    def set_file_name(self):
        """Determine filename."""
        uni_mods = os.path.join(paths.models(), 'universe/')
        self.file_name = uni_mods + 'csmd.db'

    def create_table(self):
        """Create a CSMD SQL table."""
        # Connect to database
        conn = sqlite3.connect(self.file_name)
        c = conn.cursor()

        # Set array of coordinates
        zs = np.arange(0., 6. + self.step, self.step)

        # Create database
        self.c.execute('create table csmd (z real, csmd real)')

        results = []

        H_0 = self.H_0
        W_m = self.W_m
        W_v = self.W_v

        def integral(z):
            z1 = z + 1
            return z1**1.7/(1+(z1/2.9)**5.6)*(1/(H_0*(W_m*z1**3+W_v)**0.5))

        def csmd(z):
            return 0.01095*quad(integral, z, np.inf)[0]

        vec_csmd = np.vectorize(csmd)

        # Give an update on the progress
        pprint('Creating a CSMD lookup table (only needs to be done once)')

        csmds = vec_csmd(zs)

        results = list(zip(np.around(zs, decimals=self.rounding), csmds))

        # Save results to database
        c.executemany('insert into csmd values (?,?)', results)

        # Make for easier searching
        c.execute('create index ix on csmd (z, csmd)')

        # Save
        conn.commit()

    def lookup(self, z):
        """Look up associated values with input values.

        Args:
            z (array): Redshifts

        Returns:
            array: Stellar mass density at given redshift

        """
        # Connect to database
        conn = sqlite3.connect(self.file_name)
        c = conn.cursor()

        smd = np.ones_like(z)

        # Round values
        z = np.round(z, self.rounding)

        # Search database
        query = 'select csmd from csmd where z=? limit 1'

        for i, r in enumerate(z):
            smd[i] = c.execute(query, [r]).fetchone()[0]

        # Close database
        conn.close()

        return smd
