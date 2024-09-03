"""Create a lookup tables for redshift and the NE2001, YMW16 dispersion measure."""

import os
import numpy as np
#import numexpr as ne
import time
import bisect
import sys
from scipy.integrate import quad
from tqdm import tqdm
from joblib import Parallel, delayed
import astropy.units as u
from astropy.cosmology import Planck13, Planck18, z_at_value

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
        self.rounding = 1

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
        #self.file_name = uni_mods + 'dm_mw.db'
        self.file_name = uni_mods + 'dm_mw_ne2001_1d.npy'
        if self.test:
            uni_mods = os.path.join(paths.models(), 'universe/')
            self.file_name = uni_mods + 'test_dm_mw_ne2001_1d.npy'

    def create_table(self, parallel=True):
        """Create a lookup table for dispersion measure."""
        step = 1
        gls = np.arange(-180., 180. + step, step).round(1)
        gbs = np.arange(-90., 90. + step, step).round(1)
        dist = 0.1  # [Gpc]

        gls = gls.astype(np.float32)
        gbs = gbs.astype(np.float32)

        DM_MW_Table = {}
        start = time.time()
        for gl in gls:
            for gb in gbs:
                print(gl, gb)
                if gl in DM_MW_Table:
                    DM_MW_Table[gl].update({gb: go.ne2001_dist_to_dm(dist, gl, gb)})
                else:
                    DM_MW_Table.update({gl: {gb: go.ne2001_dist_to_dm(dist, gl, gb)}})
        
        np.save(self.file_name, DM_MW_Table)
        
    def lookup(self, gal, gab):
        """Look up associated milky way dispersion measure with gal coords.

        Args:
            gl (array): Galactic longitude [fractional degrees]
            gb (array): Galactic latitude [fractional degrees]

        Returns:
            dm_mw (float): Galactic dispersion measure [pc*cm^-3]

        """
        # Load dm table
        dm_mw_table_1d = np.load(self.file_name, allow_pickle=True)

        # Round values
        gal = np.round(gal)
        gab = np.round(gab)
        index = (gal - (-180))*181 + (gab - (-90))
        #index = ne.evaluate("(gal - (-180))*181 + (gab - (-90))")
        
        dm_mw = dm_mw_table_1d[np.searchsorted(dm_mw_table_1d[:, 0], index)][:, 1]
        
        return dm_mw

class YMW16Table:
    """Create/use a NE2001 lookup table for dispersion measure."""

    def __init__(self, test=False):
        """Initializing."""
        self.test = test
        self.set_file_name()

        # Setup database
        self.db = False
        self.step = 0.1
        self.rounding = 1

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
        #self.file_name = uni_mods + 'dm_mw.db'
        self.file_name = uni_mods + 'dm_mw_ymw16_1d.npy'
        if self.test:
            uni_mods = os.path.join(paths.models(), 'universe/')
            self.file_name = uni_mods + 'test_dm_mw_ymw16_1d.npy'

    def create_table(self, parallel=True):
        """Create a lookup table for dispersion measure."""
        step = 1
        gls = np.arange(-180., 180. + step, step).round(1)
        gbs = np.arange(-90., 90. + step, step).round(1)
        dist = 0.1  # [Gpc]

        gls = gls.astype(np.float32)
        gbs = gbs.astype(np.float32)

        DM_MW_Table = {}
        start = time.time()
        for gl in gls:
            for gb in gbs:
                print(gl, gb)
                if gl in DM_MW_Table:
                    DM_MW_Table[gl].update({gb: go.ymw16_dist_to_dm(dist, gl, gb)})
                else:
                    DM_MW_Table.update({gl: {gb: go.ymw16_dist_to_dm(dist, gl, gb)}})
        
        np.save(self.file_name, DM_MW_Table)
        pprint('Finished DM table')
    
    def lookup(self, gal, gab):
        """Look up associated milky way dispersion measure with gal coords.

        Args:
            gl (array): Galactic longitude [fractional degrees]
            gb (array): Galactic latitude [fractional degrees]

        Returns:
            dm_mw (float): Galactic dispersion measure [pc*cm^-3]

        """
        # Load dm table
        dm_mw_table_1d = np.load(self.file_name, allow_pickle=True)

        gal = np.round(gal)
        gab = np.round(gab)
        index = (gal - (-180))*181 + (gab - (-90))
        #index = ne.evaluate("(gal - (-180))*181 + (gab - (-90))")
        
        dm_mw = dm_mw_table_1d[np.searchsorted(dm_mw_table_1d[:, 0], index)][:, 1]
        
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
    cdf_sfr and cdf_smd and several delayed cdf_sfr.

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

        #self.file_name = uni_mods + f'{f}.db'
        self.file_name = uni_mods + f'{f}.npy'
        
        if self.test:
            self.file_name = uni_mods + 'cosmo_test.npy'

    def create_table(self):
        """Create a lookup table for distances."""
        m = ['Creating a distance table',
             '  - Only needs to happen once',
             '  - May take up to 2m on a single core']
        for n in m:
            pprint(n)

        H_0 = self.H_0
        W_m = self.W_m
        W_v = self.W_v

        W_k = 1.0 - W_m - W_v  # Omega curvature

        if W_k != 0.0:
            pprint('Careful - Your cosmological parameters do not sum to 1.0')

        n_cpus = 128
        zs = np.arange(0, self.z_max+self.step, self.step)

        pprint('  - Calculating parameters at various redshifts')
        conv = go.Redshift(zs, H_0=H_0, W_m=W_m, W_v=W_v)
        dists = conv.dist_co()
        vols = conv.vol_co()

        # Get dV
        dvols = np.zeros_like(vols)
        dvols[1:] = np.diff(vols)

        pprint('  - Calculating Star Formation Rate')
        # Get pdf sfr
        pdf_sfr = np.array(Parallel(n_jobs=n_cpus)(delayed(sfr)(i) for i in zs))*dvols
        cdf_sfr = np.cumsum(pdf_sfr)
        cdf_sfr /= cdf_sfr[-1] # Normalize

        pprint('  - Calculating Stellar Mass Density')
        # Get pdf csmd
        pdf_smd = np.array(Parallel(n_jobs=n_cpus)(delayed(smd)(i, H_0=H_0, W_m=W_m, W_v=W_v) for i in zs))*dvols
        cdf_smd = np.cumsum(pdf_smd)
        cdf_smd /= cdf_smd[-1] # Normalize
        
        pprint('  - Calculating Delayed Star Formation Rate - 0.1 Gyr')
        # Get pdf delayed sfr 0.1 Gyr
        pdf_dsfr0d1 = np.array(Parallel(n_jobs=n_cpus)(delayed(delayed_sfr)(i, 0.1) for i in zs))*dvols
        cdf_dsfr0d1 = np.cumsum(pdf_dsfr0d1)
        cdf_dsfr0d1 /= cdf_dsfr0d1[-1] # Normalize
        
        pprint('  - Calculating Delayed Star Formation Rate - 0.5 Gyr')
        # Get pdf delayed sfr 0.5 Gyr
        pdf_dsfr0d5 = np.array(Parallel(n_jobs=n_cpus)(delayed(delayed_sfr)(i, 0.5) for i in zs))*dvols
        cdf_dsfr0d5 = np.cumsum(pdf_dsfr0d5)  # Unnormalized
        cdf_dsfr0d5 /= cdf_dsfr0d5[-1]
        
        pprint('  - Calculating Delayed Star Formation Rate - 1 Gyr')
        # Get pdf delayed sfr 1 Gyr
        pdf_dsfr1 = np.array(Parallel(n_jobs=n_cpus)(delayed(delayed_sfr)(i, 1) for i in zs))*dvols
        cdf_dsfr1 = np.cumsum(pdf_dsfr1)
        cdf_dsfr1 /= cdf_dsfr1[-1] # Normalize
        
        lookback_times = Planck18.lookback_time(zs).value

        results = np.stack((zs, dists, vols, dvols, cdf_sfr, cdf_smd, cdf_dsfr0d1, cdf_dsfr0d5, cdf_dsfr1, lookback_times)).T

        pprint('  - Saving values to database')
        
        np.save(self.file_name, results)
        
        pprint('Finished distance table')
    
    def lookup(self, z=None, dist_co=None, vol_co=None, dvol_co=None,
               cdf_sfr=None, cdf_smd=None,
               cdf_dsfr0d1=None, cdf_dsfr0d5=None, cdf_dsfr1=None,
               lookback_time=None):
        """Look up associated values with input values."""
        
        distance_table = np.load(self.file_name, allow_pickle=True)

        # Check what's being looked up, set all other keywords to same length
        kw = {'z': z,
              'dist':  dist_co,
              'vol': vol_co,
              'dvol': dvol_co,
              'cdf_sfr': cdf_sfr,
              'cdf_smd': cdf_smd,
              'cdf_dsfr0d1': cdf_dsfr0d1,
              'cdf_dsfr0d5': cdf_dsfr0d5,
              'cdf_dsfr1': cdf_dsfr1,
              'lookback_time': lookback_time}
        
        col = -1
        for key, value in kw.items():
            col += 1
            if value is not None:
                in_par = key
                break

        for key, value in kw.items():
            if key != in_par:
                kw[key] = np.ones_like(kw[in_par])

        keys = list(kw.keys())

        # Search database
        start = time.time()
        
        d = distance_table[np.searchsorted(distance_table[:, keys.index(in_par)], kw[in_par])]
        for key in keys:
            if key == in_par:
                continue
            kw[key] = d[:, keys.index(key)]

        return list(kw.values())

def delayed_sfr(z, delay_time):
    """Return the number density of star forming rate at redshift z.

    Follows Madau & Dickinson (2014), eq. 15. For more info see
    https://arxiv.org/pdf/1403.0007.pdf
    """
    #Make sure we do not exceed the z_at_value limit
    z_lim = z_at_value(Planck18.age, delay_time * u.Gyr) - 0.002
    z = np.piecewise(z, [z < z_lim, z >= z_lim], 
                        [lambda z:np.array([z_at_value(Planck18.age, 13.7869 * u.Gyr - Planck18.lookback_time(i) - delay_time * u.Gyr) for i in z]), 
                         lambda z:999])
    
    return (1+z)**2.7/(1+((1+z)/2.9)**5.6)

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
