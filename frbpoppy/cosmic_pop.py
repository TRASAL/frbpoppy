"""Class to generate a cosmic population of FRBs."""
import numpy as np

from frbpoppy.log import pprint
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
        """Set the model for generating the directions of the frb sources."""
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
            model (str): Option of 'ne2001'
        """
        if not isinstance(model, str):
            self.dm_mw_func = lambda: model(**kwargs)
            return

        # Distribution from which to draw dm_mw
        if model == 'ne2001':
            self.dm_mw_func = lambda: pc.NE2001Table().lookup(self.frbs.gl,
                                                              self.frbs.gb)
        else:
            raise ValueError('set_dm_mw input not recognised')

    def gen_dm_mw(self):
        """Generate Milky Way dispersion measure."""
        self.frbs.dm_mw = self.dm_mw_func()

    def set_dm_igm(self, model='ioka', **kwargs):
        """Set intergalactic dispersion measure model.

        Args:
            model (str): Option of 'ioka'
            slope (float): Slope of the DM-z relationship
            sigma (float): Spread around the DM-z relationship
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
            mu (float): Mean DM [pc/cm^3].
            sigma (float): Standard deviation DM [pc/cm^3].
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
            per_source (str): Model for a single source burst distribution.
                Options from ('same', 'different', 'gauss')
        If model == 'uniform':
            low (float): Minimum pulse width [ms]
            high (float): Maximum pulse width [ms]
        If model == 'lognormal':
            mu (float): Mean pulse width [ms]
            sigma (float): Standard deviation pulse width [ms]
        If per_source == 'gauss':
            src_sigma (float): Standard deviation per source [ms]

        """
        # Use your own function if you want
        if not isinstance(model, str):
            self.w_func = lambda: model(**kwargs)
            return

        # Each burst from the same source: same or different widths?
        if per_source == 'same':
            shape = self.n_srcs
        elif per_source == 'different':
            shape = self.shape
        # Or use a gaussian per source.
        elif per_source.startswith('gauss'):
            self.w_func = lambda: wd.gauss_per_source(dist=getattr(wd, model),
                                                      shape=self.shape,
                                                      z=self.frbs.z,
                                                      **kwargs)
            return

        # Distribution from which to draw pulse widths
        if model == 'uniform':
            self.w_func = lambda: wd.uniform(shape=shape, z=self.frbs.z,
                                             **kwargs)
        elif model == 'lognormal':
            self.w_func = lambda: wd.lognormal(shape=shape, z=self.frbs.z,
                                               **kwargs)
        else:
            raise ValueError('set_w input not recognised')

    def gen_w(self):
        """Generate pulse widths [ms]."""
        self.frbs.w_int, self.frbs.w_arr = self.w_func()

    def set_si(self, model='gauss', per_source='same', **kwargs):
        """Set spectral index model.

        Args:
            model (str): Options from ('gauss')
            per_source (str): Model for a single source burst distribution.
                Options from ('same', 'different', 'gauss')
        If model == 'gauss':
            mu (float): Mean spectral index
            sigma (float): Standard deviation spectral index
        If per_source == 'gauss':
            src_sigma (float): Standard deviation per source

        """
        if not isinstance(model, str):
            self.si_func = lambda: model(**kwargs)
            return

        # Each burst from the same source: same or different si?
        if per_source == 'same':
            shape = self.n_srcs
        elif per_source == 'different':
            shape = self.shape
        # Or use a gaussian per source
        elif per_source.startswith('gauss'):
            dist = getattr(sid, model)
            self.si_func = lambda: sid.gauss_per_source(dist=dist,
                                                        shape=self.shape,
                                                        **kwargs)
            return

        # Distribution from which to draw spectral indices
        if model.startswith('gauss'):
            self.si_func = lambda: sid.gauss(shape=shape, **kwargs)
        else:
            raise ValueError('set_si input not recognised')

    def gen_si(self):
        """Generate spectral indices."""
        self.frbs.si = self.si_func()

    def set_lum(self, model='powerlaw', per_source='same', **kwargs):
        """Set luminosity function [ergs/s].

        Args:
            model (str): Options from ('powerlaw')
            per_source (str): Model for a single source burst distribution.
                Options from ('same', 'different', 'gauss')

        If model == 'powerlaw':
            low (float): Minimum bolometric luminosity [ergs/s]
            high (float): Maximum bolometric luminosity [ergs/s]
            power (float): Power of luminosity function

        If per_source == 'gauss':
            src_sigma (float): Standard deviation per source

        """
        if not isinstance(model, str):
            self.lum_func = lambda: model(**kwargs)
            return

        # Each burst from the same source: same or different luminosities?
        if per_source == 'same':
            shape = self.n_srcs
        elif per_source == 'different':
            shape = self.shape
        # Or use a gaussian per source
        elif per_source.startswith('gauss'):
            dist = getattr(ld, model)
            self.lum_func = lambda: ld.gauss_per_source(dist=dist,
                                                        shape=self.shape,
                                                        **kwargs)
            return

        # Draw luminosities from a powerlaw
        if model == 'powerlaw':
            self.lum_func = lambda: ld.powerlaw(shape=shape, **kwargs)
        else:
            raise ValueError('set_lum input not recognised')

    def gen_lum(self):
        """Generate luminosities [ergs/s]."""
        self.frbs.lum_bol = self.lum_func()

    def set_time(self, model='regular', **kwargs):
        """Set model from which to generate time stamps.

        Args:
            model (str): Options from ('regular', 'clustered', 'poisson')
        If model == 'regular':
            lam (float): Number of bursts per day
        If model == 'poisson':
            lam (float): Expected number of bursts per day
        If model == 'clustered':
            r (float): Rate parameter
            k (float): Shape parameter
        """
        if not isinstance(model, str):
            # These lambda functions look complex, but aren't.
            # They merely stop the function from running immediately
            self.time_func = lambda: model(**kwargs)
            return

        if model == 'regular':
            self.time_func = lambda: td.regular(n_srcs=self.n_srcs,
                                                n_days=self.n_days,
                                                z=self.frbs.z,
                                                **kwargs)
        elif model.startswith('poisson'):
            self.time_func = lambda: td.poisson(n_srcs=self.n_srcs,
                                                n_days=self.n_days,
                                                z=self.frbs.z,
                                                **kwargs)
        elif model == 'clustered':
            self.time_func = lambda: td.clustered(n_srcs=self.n_srcs,
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
        # Generate time stamps
        self.frbs.time = self.time_func()
        # Set size for all other parameters
        self.shape = self.frbs.time.shape

    def generate(self):
        """Generate a full CosmicPopulation."""
        pprint(f'Generating {self.name} population')
        self.gen_index()
        self.gen_time()
        self.gen_dist()
        self.gen_direction()
        self.gen_gal_coords()
        self.gen_dm()
        self.gen_w()
        self.gen_lum()
        self.gen_si()
        pprint(f'Finished generating {self.name} population')

    @classmethod
    def simple(cls, n_srcs, generate=False):
        """Set up a simple, local population."""
        pop = cls(n_srcs=n_srcs, n_days=1, name='simple', repeaters=False,
                  generate=generate)
        pop.set_dist(model='vol_co', z_max=0.01, alpha=-1.5,
                     H_0=67.74, W_m=0.3089, W_v=0.6911)
        pop.set_dm(mw=False, igm=False, host=False)
        pop.set_emission_range(low=100e6, high=10e9)
        pop.set_lum(model='powerlaw', low=1e38, high=1e38, power=0)
        pop.set_w(model='uniform', low=10, high=10)
        pop.set_si(model='gauss', mu=0, sigma=0)
        return pop

    @classmethod
    def complex(cls, n_srcs, generate=False):
        """Set up a complex population."""
        pop = cls(n_srcs=n_srcs, n_days=1, name='complex', repeaters=False,
                  generate=generate)
        pop.set_dist(model='vol_co', z_max=2.5, alpha=-1.5,
                     H_0=67.74, W_m=0.3089, W_v=0.6911)
        pop.set_dm_host(model='gauss', mu=100, sigma=200)
        pop.set_dm_igm(model='ioka', slope=1000, sigma=None)
        pop.set_dm_mw(model='ne2001')
        pop.set_emission_range(low=100e6, high=10e9)
        pop.set_lum(model='powerlaw', low=1e39, high=1e45, power=0)
        pop.set_w(model='lognormal', mu=0.1, sigma=0.7)
        pop.set_si(model='gauss', mu=-1.4, sigma=1)
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
