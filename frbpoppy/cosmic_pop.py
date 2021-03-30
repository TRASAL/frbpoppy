"""Class to generate a cosmic population of FRBs."""
import numpy as np

from frbpoppy.misc import pprint
from frbpoppy.number_density import NumberDensity
from frbpoppy.population import Population
import frbpoppy.direction_dists as did
import frbpoppy.dm_dists as dmd
import frbpoppy.time_dists as td
import frbpoppy.w_dists as wd
import frbpoppy.si_dists as sid
import frbpoppy.lum_dists as ld
import frbpoppy.galacticops as go
import frbpoppy.precalc as pc
import frbpoppy.healpix_maps as hpm

class CosmicPopulation(Population):
    """Generate a cosmic FRB population."""

    def __init__(self,
                 n_srcs=1e4,
                 n_days=1,
                 name='cosmic',
                 repeaters=False,
                 generate=False):
        """Generate a popuation of FRBs.

        Args:
            n_srcs (int): Number of FRB sources to generate.
            n_days (float): Number of days over which FRBs are to be generated.
            name (str): Population name.
            repeaters (bool): Whether to generate a repeater population.
            generate (bool): Whether to create a population.

        Returns:
            Population: Population of FRBs

        """
        # Set up general population arguments
        Population.__init__(self)

        # Give population a name
        self.name = 'cosmic'
        if name:
            self.name = name

        # Set population arguments
        self.n_srcs = int(n_srcs)
        self.n_days = n_days
        self.repeaters = repeaters
        self.shape = (self.n_srcs,)

        # If wanting repeaters
        if self.repeaters:
            self.set_time()

        # Set up default models
        self.set_emission_range()
        self.set_dist()
        self.set_direction()
        self.set_lum()
        self.set_si()
        self.set_w()
        self.set_dm_mw()
        self.set_dm_igm()
        self.set_dm_host()
        self.set_dm()

        # Whether to start generating a Cosmic Population
        if generate:
            self.generate()

    def gen_index(self):
        """Generate indices for over each FRB source."""
        self.frbs.index = np.arange(self.n_srcs)

    def set_emission_range(self, low=100e6, high=10e9):
        """Set the emission range [Hz].

        The frequency range between which FRB sources should emit the given
        bolometric luminosity.

        Args:
            f_min (float): Lowest source frequency [Hz].
            f_max (float): Highest source frequency [Hz].
        """
        self.f_min = low
        self.f_max = high

    def gen_precalc(self):
        """Check whether pre-calculations have been run."""
        pc.DistanceTable(H_0=self.H_0, W_m=self.W_m, W_v=self.W_v)

    def set_dist(self, model='vol_co', **kwargs):
        """Set the number density model for calculating distances.

        Args:
            model (str): Number density model to use. Choice from
                ('vol_co', 'sfr', 'smd').
            z_max (float): Maximum redshift.
            H_0 (float): Hubble constant.
            W_m (float): Density parameter Ω_m.
            W_v (float): Cosmological constant Ω_Λ.
            alpha (float): Desired log N log S slope for a perfect,
                non-cosmological population.
        """
        # Option to use your own model
        if not isinstance(model, str):
            self.dist_func = lambda: model(**kwargs)
            return

        # I sometimes use 'constant' instead of 'vol_co'
        if model == 'constant':
            model = 'vol_co'

        # Check whether recognised number density model
        if model not in ['vol_co', 'sfr', 'smd']:
            raise ValueError('set_dist input not recognised')

        # Set number density model
        # Don't fear the lambda, merely delays executing the function
        self.dist_func = lambda: NumberDensity(model=model, **kwargs)

    def gen_dist(self):
        """Generate source distances."""
        n_model = self.dist_func()
        self.vol_co_max = n_model.vol_co_max
        self.frbs.z, self.frbs.dist_co = n_model.draw(self.n_srcs)

    def gen_gal_coords(self):
        """Generate galactic coordinates."""
        frbs = self.frbs
        # Get the proper distance
        dist_pr = frbs.dist_co/(1+frbs.z)
        # Convert into galactic coordinates
        frbs.gx, frbs.gy, frbs.gz = go.lb_to_xyz(frbs.gl, frbs.gb, dist_pr)

    def set_direction(self, model='uniform', **kwargs):
        """Set the model for generating the directions of the frb sources.

        Args:
            model (str): Choice from ('uniform').
        if model == 'uniform':
            min_ra (float): Minimum right ascenion [frac deg].
            max_ra (float): Maximum right ascenion [frac deg].
            min_dec (float): Minimum declination [frac deg].
            max_dec (float): Maximum declination [frac deg].
        """
        # Use your own function
        if not isinstance(model, str):
            self.direction_func = lambda: model(**kwargs)
            return

        # Or use a uniform distribution
        if model == 'uniform':
            self.direction_func = lambda: did.uniform(n_srcs=self.n_srcs,
                                                      **kwargs)
        else:
            raise ValueError('set_direction input not recognised')

    def gen_direction(self):
        """Generate the direction of frbs."""
        frbs = self.frbs
        # Calculate right ascenion and declination
        frbs.ra, frbs.dec = self.direction_func()
        # Convert to galactic lat/long coordinates
        frbs.gl, frbs.gb = go.radec_to_lb(frbs.ra, frbs.dec, frac=True)

    def set_dm_mw(self, model='ne2001', **kwargs):
        """Set the model for the Milky Way dispersion measure.

        Args:
            model (str): Option of 'ne2001' or 'ymw16'.
        """
        if not isinstance(model, str):
            self.dm_mw_func = lambda: model(**kwargs)
            return

        # Distribution from which to draw dm_mw
        if model in('ne2001', 'ymw16'):
            self.dm_mw_func = lambda: hpm.dm_mw(self.frbs.gl, self.frbs.gb, model)
        else:
            raise ValueError('set_dm_mw input not recognised')

    def gen_dm_mw(self):
        """Generate Milky Way dispersion measure."""
        self.frbs.dm_mw = self.dm_mw_func()

    def set_dm_igm(self, model='ioka', **kwargs):
        """Set intergalactic dispersion measure model.

        Args:
            model (str): Option of 'ioka'.
        if model == 'ioka':
            slope (float): Slope of the DM-z relationship.
            std (float): Spread around the DM-z relationship.
            spread_dist (str): 'normal' or 'lognormal'.
        """
        # Possibility to use your own function
        if not isinstance(model, str):
            self.dm_igm_func = lambda: model(**kwargs)
            return

        # Distribution from which to draw intergalactic dm
        if model == 'ioka':
            self.dm_igm_func = lambda: dmd.ioka(z=self.frbs.z, **kwargs)
        else:
            raise ValueError('set_dm_igm input not recognised')

    def gen_dm_igm(self):
        """Generate intergalactic dispersion measure."""
        self.frbs.dm_igm = self.dm_igm_func()

    def set_dm_host(self, model='gauss', **kwargs):
        """Set host galaxy dispersion measure.

        Args:
            model (str): Options from ('gauss', 'lognormal').

        if model in ('gauss', 'lognormal'):
            mean (float): Mean DM [pc/cm^3].
            std (float): Standard deviation DM [pc/cm^3].
        if model == 'constant':
            value (float): Value to adopt [pc/cm^3].
        """
        if not isinstance(model, str):
            self.dm_host_func = lambda: model(**kwargs)
            return

        # Distribution from which to draw host dispersion measure
        if model.startswith('gauss'):
            self.dm_host_func = lambda: dmd.gauss(z=self.frbs.z,
                                                  n_srcs=self.n_srcs,
                                                  **kwargs)
        elif model == 'lognormal':
            self.dm_host_func = lambda: dmd.lognormal(z=self.frbs.z,
                                                      n_srcs=self.n_srcs,
                                                      **kwargs)
        elif model == 'constant':
            self.dm_host_func = lambda: dmd.constant(n_srcs=self.n_srcs,
                                                     **kwargs)
        else:
            raise ValueError('set_dm_host input not recognised')

    def gen_dm_host(self):
        """Generate host dispersion measure."""
        self.frbs.dm_host = self.dm_host_func()

    def set_dm(self, mw=True, igm=True, host=True):
        """Set total dispersion measure.

        Args:
            mw (bool): Whether to include a Milky Way component
            igm (bool): Whether to include an IGM component
            host (bool): Whether to include a host galaxy component
        """
        # Which components to include
        self.dm_components = []
        if mw:
            self.dm_components.append(self.gen_dm_mw)
        if igm:
            self.dm_components.append(self.gen_dm_igm)
        if host:
            self.dm_components.append(self.gen_dm_host)

        # Save those components to execute at a later stage
        def run_dm():
            [c() for c in self.dm_components]
            return self.frbs.dm_mw + self.frbs.dm_igm + self.frbs.dm_host

        self.dm_func = run_dm

    def gen_dm(self):
        """Generate total dispersion measure."""
        self.frbs.dm = self.dm_func()

    def set_w(self, model='uniform', per_source='same', **kwargs):
        """Set intrinsic pulse widths model [ms].

        Args:
            model (str): Options from ('uniform', 'lognormal')
            per_source (str): Model for a single source burst
                distribution. Options from 'same' or 'different'
        If model == 'constant':
            value (float): Pulse width [ms].
        If model == 'uniform':
            low (float): Minimum pulse width [ms].
            high (float): Maximum pulse width [ms].
        If model == 'lognormal':
            mean (float): Mean pulse width [ms].
            std (float): Standard deviation pulse width [ms].

        """
        # Each burst from the same source: same or different widths?
        if per_source == 'same':
            self.w_shape = lambda: self.n_srcs
        elif per_source == 'different':
            self.w_shape = lambda: self.shape

        # Distribution from which to draw pulse widths

        # Find available distributions to draw from
        funcs = [d for d in dir(wd) if hasattr(getattr(wd, d), '__call__')]
        funcs.remove('calc_w_arr')

        # Set function
        if model in funcs:
            func = getattr(wd, model)

            # If you're getting fancy with combined distributions
            # See examples/adapting_population_parameters.py
            self._transpose_w = False
            for kw_value in kwargs.values():
                if isinstance(kw_value, (list, np.ndarray)):
                    self.w_shape = lambda: self.shape[::-1]
                    self._transpose_w = True

            self.w_func = lambda x: func(shape=x, z=self.frbs.z, **kwargs)
        else:
            raise ValueError('set_w input model not recognised')

    def gen_w(self):
        """Generate pulse widths [ms]."""
        shape = self.w_shape()
        self.frbs.w_int, self.frbs.w_arr = self.w_func(shape)

        # From combined distribution inputs
        if self._transpose_w:
            self.frbs.w_int = self.frbs.w_int.T
            self.frbs.w_arr = self.frbs.w_arr.T

    def set_si(self, model='gauss', per_source='same', **kwargs):
        """Set spectral index model.

        Args:
            model (str): Options from ('gauss')
            per_source (str): Model for a single source burst
                distribution. Options from ('same', 'different')
        If model == 'constant':
            value (float): Default spectal index.
        If model == 'gauss':
            mean (float): Mean spectral index
            std (float): Standard deviation spectral index

        """
        # Each burst from the same source: same or different si?
        if per_source == 'same':
            self.si_shape = lambda: self.n_srcs
        elif per_source == 'different':
            self.si_shape = lambda: self.shape

        # Find available distributions to draw from
        funcs = [d for d in dir(sid) if hasattr(getattr(sid, d), '__call__')]

        # Set function
        if model in funcs:
            func = getattr(sid, model)

            # If you're getting fancy with combined distributions
            self._transpose_si = False
            for kw_value in kwargs.values():
                if isinstance(kw_value, (list, np.ndarray)):
                    self.si_shape = lambda: self.shape[::-1]
                    self._transpose_si = True

            # Distribution from which to draw spectral indices
            self.si_func = lambda x: func(shape=x, **kwargs)
        else:
            raise ValueError('set_si input not recognised')

    def gen_si(self):
        """Generate spectral indices."""
        shape = self.si_shape()
        self.frbs.si = self.si_func(shape)

        if self._transpose_si:
            self.frbs.si = self.frbs.si.T

    def set_lum(self, model='powerlaw', per_source='same', **kwargs):
        """Set luminosity function [ergs/s].

        Args:
            model (str): Options from ('powerlaw')
            per_source (str): Model for a single source burst
                distribution. Options from ('same', 'different')

        If model == 'powerlaw':
            low (float): Minimum bolometric luminosity [ergs/s]
            high (float): Maximum bolometric luminosity [ergs/s]
            power (float): Power of luminosity function

        If model == 'constant':
            value (float): Value for standard candle [ergs/s]

        """
        # Each burst from the same source: same or different luminosities?
        if per_source == 'same':
            self.lum_shape = lambda: self.n_srcs
        elif per_source == 'different':
            self.lum_shape = lambda: self.shape

        # Find available distributions to draw from
        funcs = [d for d in dir(ld) if hasattr(getattr(ld, d), '__call__')]

        # Set function
        if model in funcs:
            func = getattr(ld, model)

            # Help out the user
            for s in ['slope', 'index']:
                if s in kwargs:
                    kwargs['power'] = kwargs.pop(s)

            # If you're getting fancy with combined distributions
            self._transpose_lum = False
            for kw_value in kwargs.values():
                if isinstance(kw_value, (list, np.ndarray)):
                    self.lum_shape = lambda: self.shape[::-1]
                    self._transpose_lum = True

            # Distribution from which to draw luminosities
            self.lum_func = lambda x: func(shape=x, **kwargs)
        else:
            raise ValueError('set_lum input not recognised')

    def gen_lum(self):
        """Generate luminosities [ergs/s]."""
        shape = self.lum_shape()
        self.frbs.lum_bol = self.lum_func(shape)

        # You need multiple luminosities if repeaters
        if self.repeaters and self.frbs.lum_bol.ndim == 1:
            repeat_lums = [self.frbs.lum_bol[..., np.newaxis]]*self.shape[1]
            self.frbs.lum_bol = np.concatenate(repeat_lums, axis=1)

        if self._transpose_lum:
            self.frbs.lum_bol = self.frbs.lum_bol.T

    def set_time(self, model='regular', **kwargs):
        """Set model from which to generate time stamps.

        Args:
            model (str): Options from ('single', 'regular', 'clustered',
                'poisson', 'cyclic')
        If model == 'regular':
            rate (float): Number of bursts per day
        If model == 'poisson':
            rate (float): Expected number of bursts per day
        If model == 'clustered':
            r (float): Rate parameter
            k (float): Shape parameter
        If model == 'cyclic':
            rate (float): Number of bursts per day
            period (float): Period of activity cycle (days)
            frac (float): Fraction of activity cycle a source is active
        """
        if not isinstance(model, str):
            # These lambda functions look complex, but aren't.
            # They merely stop the function from running immediately
            self.time_func = lambda: model(**kwargs)
            return

        # Find available distributions
        funcs = [d for d in dir(td) if hasattr(getattr(td, d), '__call__')]
        internal = ['gamma', 'iteratively_gen_times', '_weibull_dist',
                    '_poisson_dist']
        for f in internal:
            funcs.remove(f)

        # Set function
        if model in funcs:
            func = getattr(td, model)

            # Distribution from which to draw time stamps
            self.time_func = lambda: func(n_srcs=self.n_srcs,
                                          n_days=self.n_days,
                                          z=self.frbs.z,
                                          **kwargs)
        else:
            raise ValueError('set_time input not recognised')

    def gen_time(self):
        """Generate time stamps."""
        # Only relevant for repeaters
        if not self.repeaters:
            return

        pprint('Adding burst times')
        self.frbs.time = self.time_func()
        # Set size for all other parameters
        self.shape = self.frbs.time.shape
        pprint('Finished adding burst times')

    def generate(self):
        """Generate a full CosmicPopulation."""
        pprint(f'Generating {self.name} population')
        self.gen_index()
        self.gen_dist()
        self.gen_time()
        self.gen_direction()
        self.gen_gal_coords()
        self.gen_dm()
        self.gen_w()
        self.gen_lum()
        self.gen_si()
        pprint(f'Finished generating {self.name} population')

    @classmethod
    def simple(cls, n_srcs, n_days=1, repeaters=False, generate=False):
        """Set up a simple, local population."""
        pop = cls(n_srcs=n_srcs, n_days=n_days, name='simple',
                  repeaters=repeaters, generate=False)
        pop.set_dist(model='vol_co', z_max=0.01, alpha=-1.5,
                     H_0=67.74, W_m=0.3089, W_v=0.6911)
        pop.set_dm(mw=False, igm=False, host=False)
        pop.set_emission_range(low=10e7, high=10e9)
        pop.set_lum(model='constant', value=1e38)
        pop.set_w(model='constant', value=10)
        pop.set_si(model='constant', value=0)
        if pop.repeaters:
            pop.set_time(model='regular', rate=2)
        if generate:
            pop.generate()
        return pop

    @classmethod
    def complex(cls, n_srcs, n_days=1, repeaters=False, generate=False):
        """Set up a complex population."""
        pop = cls(n_srcs=n_srcs, n_days=n_days, name='complex',
                  repeaters=repeaters, generate=False)
        pop.set_dist(model='vol_co', z_max=1, alpha=-1.5,
                     H_0=67.74, W_m=0.3089, W_v=0.6911)
        pop.set_dm_host(model='gauss', mean=100, std=200)
        pop.set_dm_igm(model='ioka', slope=1000, std=None)
        pop.set_dm_mw(model='ne2001')
        pop.set_dm(mw=True, igm=True, host=True)
        pop.set_emission_range(low=10e7, high=10e9)
        pop.set_lum(model='powerlaw', low=1e40, high=1e45, power=0)
        pop.set_w(model='lognormal', mean=0.1, std=1)
        pop.set_si(model='gauss', mean=-1.4, std=1)
        if pop.repeaters:
            pop.set_time(model='poisson', rate=9)
        if generate:
            pop.generate()
        return pop

    @classmethod
    def optimal(cls, n_srcs, n_days=1, repeaters=False, generate=False):
        """Set up a complex population."""
        pop = cls(n_srcs=n_srcs, n_days=n_days, name='optimal',
                  repeaters=repeaters, generate=False)
        pop.set_dist(model='vol_co', z_max=2.5, alpha=-2.2,
                     H_0=67.74, W_m=0.3089, W_v=0.6911)
        pop.set_dm_host(model='constant', value=50)
        pop.set_dm_igm(model='ioka', slope=1000, std=None)
        pop.set_dm_mw(model='ne2001')
        pop.set_dm(mw=True, igm=True, host=True)
        pop.set_emission_range(low=10e7, high=10e9)
        pop.set_lum(model='powerlaw', low=1e40, high=1e45, power=-0.8)
        pop.set_w(model='lognormal', mean=6.3e-3, std=.6)
        pop.set_si(model='constant', value=-0.4)
        if pop.repeaters:
            pop.set_time(model='poisson', rate=9)
        if generate:
            pop.generate()
        return pop

if __name__ == '__main__':
    # Quick test whether everything seems to be working or not
    import os
    import matplotlib.pyplot as plt
    pop = CosmicPopulation(1e4)
    pop.generate()

    frbs = pop.frbs

    for arg in frbs.__dict__:
        pprint(f'Plotting {arg}')
        values = getattr(frbs, arg)
        if values is not None:
            plt.hist(values, bins=50)
            plt.xlabel(arg)
            p = f'../tests/plots/{arg}.png'
            p = os.path.join(os.path.abspath(os.path.dirname(__file__)), p)
            plt.savefig(p)
            plt.clf()
