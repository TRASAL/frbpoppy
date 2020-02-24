"""Create a lookup tables for redshift and the NE2001 dispersion measure."""

import os
import numpy as np
import sqlite3
import sys
from scipy.integrate import quad
from tqdm import tqdm
from joblib import Parallel, delayed

import frbpoppy.galacticops as go
from frbpoppy.misc import pprint
from frbpoppy.paths import paths


class NE2001Table:
    """Create/use a NE2001 lookup table for dispersion measure."""

    def __init__(self, test=False):
        """Initializing."""
        self.test = test
        self.set_file_name()

        # Setup database
        self.db = False
        self.step = 0.1
        self.rounding = 2

        # For parallel processes
        self.temp_path = None

        if self.test:
            self.step = 0.1
            if os.path.exists(self.file_name):
                os.remove(self.file_name)

        if os.path.exists(self.file_name) and self.test is False:
            self.db = True
        else:
            # Calculations take quite some time
            # Provide a way for people to quit
            try:
                self.create_table()
            except KeyboardInterrupt:
                pprint('Losing all progress in calculations')
                os.remove(self.file_name)
                if self.temp:
                    os.remove(self.temp_path)
                sys.exit()

    def set_file_name(self):
        """Determine filename."""
        uni_mods = os.path.join(paths.models(), 'universe/')
        self.file_name = uni_mods + 'dm_mw.db'

        if self.test:
            uni_mods = os.path.join(paths.models(), 'universe/')
            self.file_name = uni_mods + 'test_dm_mw.db'

    def create_table(self, parallel=True):
        """Create a lookup table for dispersion measure."""
        # Connect to database
        conn = sqlite3.connect(self.file_name)
        c = conn.cursor()

        # Set array of coordinates
        gls = np.arange(-180., 180. + self.step, self.step).round(1)
        gbs = np.arange(-90., 90. + self.step, self.step).round(1)
        dist = 0.1  # [Gpc]

        gls = gls.astype(np.float32)
        gbs = gbs.astype(np.float32)

        # Create database
        c.execute('create table dm ' +
                  '(gl real, gb real, dm_mw real)')

        # Give an update on the progress
        m = ['Creating a DM lookup table',
             '  - Only needs to happen once',
             '  - Unfortunately pretty slow',
             '  - Prepare to wait for ~1.5h (4 cores)',
             '  - Time given as [time_spent<time_left] in (hh:)mm:ss',
             'Starting to calculate DM values']
        for n in m:
            pprint(n)

        n_opt = len(gls)*len(gbs)
        options = np.array(np.meshgrid(gls, gbs)).T.reshape(-1, 2)
        dm_mw = np.zeros(len(options)).astype(np.float32)

        def dm_tot(i, dm_mw):
            gl, gb = options[i]
            dm_mw[i] = go.ne2001_dist_to_dm(dist, gl, gb)

        if parallel:

            temp_path = os.path.join(paths.models(), 'universe/') + 'temp.mmap'
            self.temp_path = temp_path

            # Make a temp memmap to have a sharedable memory object
            temp = np.memmap(temp_path, dtype=dm_mw.dtype,
                             shape=len(dm_mw),
                             mode='w+')

            # Parallel process in order to populate array
            r = range(n_opt)
            j = min([4, os.cpu_count() - 1])
            print(os.cpu_count())
            Parallel(n_jobs=j)(delayed(dm_tot)(i, temp) for i in tqdm(r))

            # Map results
            r = np.concatenate((options, temp[:, np.newaxis]), axis=1)
            results = map(tuple, r.tolist())

            # Delete the temporary directory and contents
            try:
                os.remove(temp_path)
            except FileNotFoundError:
                print(f'Unable to remove {temp_path}')

        else:
            for i in tqdm(range(n_opt)):
                dm_tot(i, dm_mw)

            # Save results to database
            dm_mw = dm_mw.astype(np.float32)
            r = np.concatenate((options, dm_mw[:, np.newaxis]), axis=1)
            results = map(tuple, r.tolist())

        pprint('  - Saving results')
        c.executemany('insert into dm values (?,?,?)', results)

        # Make for easier searching
        c.execute('create index ix on dm (gl, gb)')

        # Save
        conn.commit()

        pprint('Finished DM table')

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
            dm_mw[i] = c.execute(query, [str(gl), str(gab[i])]).fetchone()[0]

        # Close database
        conn.close()

        return dm_mw


class DistanceTable:
    """
    Create/use a lookup table for comoving distance, volume, redshift etc.

    Create a list of tuples to lookup the corresponding redshift for a comoving
    distance [Gpc] (or the other way around). Uses formulas from
    Hoggs et al. (1999) for the cosmological calculations. To avoid long
    calculation times, it will check if a previous run with the same parameters
    has been done, which it will then load it. If not, it will calculate a new
    table, and save the table for later runs. Covers z, dist, vol, dvol,
    cdf_sfr and cdf_smd.

    Args:
        H_0 (float, optional): Hubble parameter. Defaults to 67.74 km/s/Mpc
        W_m (float, optional): Omega matter. Defaults to 0.3089
        W_k (float, optional): Omega vacuum. Defaults to 0.6911

    """

    def __init__(self, H_0=67.74, W_m=0.3089, W_v=0.6911, test=False):
        """Initializing."""
        self.H_0 = H_0
        self.W_m = W_m
        self.W_v = W_v
        self.test = test

        self.set_file_name()

        # Setup database
        self.db = False
        self.step = 0.00001
        self.z_max = 6.5

        if self.test:
            self.step = 0.001
            self.z_max = 6.5
            if os.path.exists(self.file_name):
                os.remove(self.file_name)

        if os.path.exists(self.file_name) and self.test is False:
            self.db = True
        else:
            # Calculations take quite some time
            # Provide a way for people to quit
            try:
                self.create_table()
            except KeyboardInterrupt:
                pprint('Losing all progress in calculations')
                os.remove(self.file_name)
                sys.exit()

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

        self.file_name = uni_mods + f'{f}.db'

        if self.test:
            self.file_name = uni_mods + 'cosmo_test.db'

    def create_table(self):
        """Create a lookup table for distances."""
        m = ['Creating a distance table',
             '  - Only needs to happen once',
             '  - May take up to 2m on a single core']
        for n in m:
            pprint(n)

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
        t = 'real'
        par = f'(z {t}, dist {t}, vol {t}, dvol {t}, cdf_sfr {t}, cdf_smd {t})'
        s = f'create table distances {par}'
        c.execute(s)

        results = []

        pprint('  - Calculating parameters at various redshifts')
        conv = go.Redshift(zs, H_0=H_0, W_m=W_m, W_v=W_v)
        dists = conv.dist_co()
        vols = conv.vol_co()

        # Get dV
        dvols = np.zeros_like(vols)
        dvols[1:] = np.diff(vols)

        pprint('  - Calculating Star Formation Rate')
        # Get pdf sfr
        pdf_sfr = sfr(zs)*dvols
        cdf_sfr = np.cumsum(pdf_sfr)  # Unnormalized
        cdf_sfr /= cdf_sfr[-1]

        pprint('  - Calculating Stellar Mass Density')
        # Get pdf csmd
        pdf_smd = smd(zs, H_0=H_0, W_m=W_m, W_v=W_v)*dvols
        cdf_smd = np.cumsum(pdf_smd)  # Unnormalized
        cdf_smd /= cdf_smd[-1]

        results = np.stack((zs, dists, vols, dvols, cdf_sfr, cdf_smd)).T

        pprint('  - Saving values to database')
        # Save results to database
        data = map(tuple, results.tolist())
        c.executemany('insert into distances values (?,?,?,?,?,?)', data)

        # Make for easier searching
        # I don't really understand SQL index names...
        c.execute('create index ix on distances (z)')
        c.execute('create index ixx on distances (dist)')
        c.execute('create index ixxx on distances (vol)')
        c.execute('create index ixxxx on distances (dvol)')
        c.execute('create index ixxxxx on distances (cdf_sfr)')
        c.execute('create index ixxxxxx on distances (cdf_smd)')

        # Save
        conn.commit()

        pprint('Finished distance table')

    def lookup(self, z=None, dist_co=None, vol_co=None, dvol_co=None,
               cdf_sfr=None, cdf_smd=None):
        """Look up associated values with input values."""
        # Connect to database
        conn = sqlite3.connect(self.file_name)
        c = conn.cursor()

        # Check what's being looked up, set all other keywords to same length
        kw = {'z': z,
              'dist':  dist_co,
              'vol': vol_co,
              'dvol': dvol_co,
              'cdf_sfr': cdf_sfr,
              'cdf_smd': cdf_smd}

        for key, value in kw.items():
            if value is not None:
                in_par = key
                break

        for key, value in kw.items():
            if key != in_par:
                kw[key] = np.ones_like(kw[in_par])

        keys = list(kw.keys())

        # Search database
        query = f'select * from distances where {in_par} > ? limit 1'

        for i, r in enumerate(kw[in_par]):
            d = c.execute(query, [str(r)]).fetchone()
            for ii, key in enumerate(keys):
                if key == in_par:
                    continue

                kw[key][i] = d[ii]

        # Close database
        conn.close()

        return list(kw.values())


def sfr(z):
    """Return the number density of star forming rate at redshift z.

    Follows Madau & Dickinson (2014), eq. 15. For more info see
    https://arxiv.org/pdf/1403.0007.pdf
    """
    return (1+z)**2.7/(1+((1+z)/2.9)**5.6)


def smd(z, H_0=67.74, W_m=0.3089, W_v=0.6911):
    """Return the number density of Stellar Mass Density at redshift z.

    Follows Madau & Dickinson (2014), eq. 2 & 15. For more info see
    https://arxiv.org/pdf/1403.0007.pdf
    """
    def integral(z):
        z1 = z + 1
        return z1**1.7/(1+(z1/2.9)**5.6)*(1/(H_0*(W_m*z1**3+W_v)**0.5))

    def csmd(z):
        return 0.01095*quad(integral, z, np.inf)[0]

    vec_csmd = np.vectorize(csmd)

    return vec_csmd(z)
