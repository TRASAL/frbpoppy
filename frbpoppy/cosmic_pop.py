"""Class to generate a cosmic population of FRBs."""
import numpy as np

from frbpoppy.log import pprint
from frbpoppy.number_density import NumberDensity
from frbpoppy.population import Population
import frbpoppy.distributions as dis
import frbpoppy.galacticops as go
import frbpoppy.precalc as pc


class CosmicPopulation(Population):
    """Generate a cosmic FRB population."""

    def __init__(self,
                 n_gen,
                 days=1,
                 name='cosmic',
                 H_0=67.74,
                 W_m=0.3089,
                 W_v=0.6911,
                 dm_host_model='normal',
                 dm_host_mu=100,
                 dm_host_sigma=200,
                 dm_igm_index=1000,
                 dm_igm_sigma=None,
                 dm_mw_model='ne2001',
                 emission_range=[10e6, 10e9],
                 lum_range=[1e40, 1e45],
                 lum_index=0,
                 n_model='sfr',
                 alpha=-1.5,
                 w_model='lognormal',
                 w_range=[0.1, 10],
                 w_mu=0.1,
                 w_sigma=0.5,
                 si_mu=-1.4,
                 si_sigma=1.,
                 z_max=2.5,
                 generate=True):
        """Generate a popuation of FRBs.

        Args:
            n_gen (int): Number of FRB sources/sky/time to generate.
            days (float): Number of days over which FRBs are generated.
            name (str): Population name.
            H_0 (float): Hubble constant.
            W_m (float): Density parameter Ω_m.
            W_v (float): Cosmological constant Ω_Λ.
            dm_host_model (float): Dispersion measure host model. Options are
                'normal' or 'lognormal'.
            dm_host_mu (float): Mean dispersion measure host [pc/cm^3].
            dm_host_sigma (float): Deviation dispersion measure host [pc/cm^3].
            dm_igm_index (float): Dispersion measure slope for IGM [pc/cm^3].
            dm_igm_sigma (float): Scatter around dm_igm. Defaults 0.2*slope*z
            dm_mw_model (str): Dispersion measure model for the Milky Way.
                Options are 'ne2001' or 'zero'.
            emission_range (list): The frequency range [Hz] between which FRB
                sources should emit the given bolometric luminosity.
            lum_range (list): Bolometric luminosity (distance) range [erg/s].
            lum_index (float): Power law index.
            n_model (str): Number density model. Either 'vol_co', 'sfr' or
                'smd'.
            alpha (float): Desired logN/logS of perfectly detected population.
            w_model (str): Pulse width model, 'lognormal' or 'uniform'.
            w_range (list): Pulse width range [ms].
            w_mu (float): Mean pulse width [ms].
            w_sigma (float): Deviation pulse width [ms].
            si_mu (float): Mean spectral index.
            si_sigma (float): Standard deviation spectral index.
            z_max (float): Maximum redshift.
            generate (bool): Whether to create a population

        Returns:
            Population: Population of FRBs.

        """
        # Set up population
        Population.__init__(self)
        self.alpha = alpha
        self.dm_host_model = dm_host_model
        self.dm_host_mu = dm_host_mu
        self.dm_host_sigma = dm_host_sigma
        self.dm_igm_index = dm_igm_index
        self.dm_igm_sigma = dm_igm_sigma
        self.dm_mw_model = dm_mw_model
        self.f_max = emission_range[1]
        self.f_min = emission_range[0]
        self.H_0 = H_0
        self.lum_max = lum_range[1]
        self.lum_min = lum_range[0]
        self.lum_pow = lum_index
        self.name = name
        self.n_gen = int(n_gen)
        self.n_model = n_model
        self.si_mu = si_mu
        self.si_sigma = si_sigma
        self.time = days * 86400  # Convert to seconds
        self.w_model = w_model
        self.w_max = w_range[1]
        self.w_min = w_range[0]
        self.w_mu = w_mu
        self.w_sigma = w_sigma
        self.W_m = W_m
        self.W_v = W_v
        self.z_max = z_max

        # Whether to start generating a Cosmic Population
        if generate:
            self.generate()

    def gen_dist(self):
        """Generate distances."""
        # Cosmology calculations
        r = go.Redshift(self.z_max,
                        H_0=self.H_0,
                        W_m=self.W_m,
                        W_v=self.W_v)
        self.dist_co_max = r.dist_co()
        self.vol_co_max = r.vol_co()

        # Ensure precalculations are done if necessary
        pc.DistanceTable(H_0=self.H_0, W_m=self.W_m, W_v=self.W_v)

        # Set up number density
        n_den = NumberDensity(model=self.n_model,
                              z_max=self.z_max,
                              alpha=self.alpha,
                              H_0=self.H_0,
                              W_m=self.W_m,
                              W_v=self.W_v).draw

        frbs = self.frbs

        # Draw from number density
        frbs.z, frbs.dist_co = n_den(self.n_gen)

    def gen_direction(self):
        """Generate the direction of frbs."""
        frbs = self.frbs

        # Add random directional coordinates
        u = np.random.uniform
        frbs.ra = u(0, 360, self.n_gen)
        frbs.dec = np.rad2deg(np.arccos(u(-1, 1, self.n_gen))) - 90

        # Convert to galactic coordinates
        frbs.gl, frbs.gb = go.radec_to_lb(frbs.ra, frbs.dec, frac=True)

    def gen_gal_coords(self):
        """Generate galactic coordinates."""
        frbs = self.frbs
        # Get the proper distance
        dist_pr = frbs.dist_co/(1+frbs.z)

        # Convert into galactic coordinates
        frbs.gx, frbs.gy, frbs.gz = go.lb_to_xyz(frbs.gl, frbs.gb, dist_pr)

    def gen_dm_host(self):
        """Generate dm host contributions."""
        frbs = self.frbs
        # Dispersion measure of the host (Tendulkar)
        if self.dm_host_model == 'normal':
            frbs.dm_host = dis.trunc_norm(self.dm_host_mu,
                                          self.dm_host_sigma,
                                          self.n_gen).astype(np.float32)
        elif self.dm_host_model == 'lognormal':
            frbs.dm_host = np.random.lognormal(self.dm_host_mu,
                                               self.dm_host_sigma,
                                               self.n_gen).astype(np.float32)

        frbs.dm_host = frbs.dm_host / (1 + frbs.z)

    def gen_dm(self):
        """Generate dispersion measures."""
        frbs = self.frbs

        # Dispersion measure of the Milky Way
        if self.dm_mw_model == 'ne2001':
            frbs.dm_mw = pc.NE2001Table().lookup(frbs.gl, frbs.gb)
        elif self.dm_mw_model == 'zero':
            frbs.dm_mw = np.zeros_like(frbs.z)

        # Dispersion measure of the intergalactic medium
        frbs.dm_igm = go.ioka_dm_igm(frbs.z,
                                     slope=self.dm_igm_index,
                                     sigma=self.dm_igm_sigma)

        # Dispersion measure of the host (Tendulkar)
        self.gen_dm_host()

        # Total dispersion measure
        frbs.dm = frbs.dm_mw + frbs.dm_igm + frbs.dm_host

    def gen_w(self, shape):
        """Generate pulse widths."""
        frbs = self.frbs
        # Get a random intrinsic pulse width [ms]
        if self.w_model == 'lognormal':
            frbs.w_int = np.random.lognormal(self.w_mu, self.w_sigma,
                                             shape).astype(np.float32)

        if self.w_model == 'uniform':
            frbs.w_int = np.random.uniform(self.w_min, self.w_max,
                                           shape).astype(np.float32)

        # Calculate the pulse width upon arrival to Earth
        if isinstance(shape, tuple):
            frbs.w_arr = frbs.w_int*(1+frbs.z[:, None])
        else:
            frbs.w_arr = frbs.w_int*(1+frbs.z)

    def gen_lum(self, shape):
        """Generate luminosities."""
        frbs = self.frbs
        # Add bolometric luminosity [erg/s]
        frbs.lum_bol = dis.powerlaw(self.lum_min,
                                    self.lum_max,
                                    self.lum_pow,
                                    shape).astype(np.float64)

    def gen_si(self, shape):
        """Generate spectral indices."""
        frbs = self.frbs
        # Add spectral index
        frbs.si = np.random.normal(self.si_mu, self.si_sigma,
                                   shape).astype(np.float32)

    def generate(self):
        """Generate all manner of intrinsic parameters."""
        # Let user know what's happening
        pprint(f'Generating {self.name} population')
        self.gen_dist()
        self.gen_direction()
        self.gen_gal_coords()
        self.gen_dm()
        self.gen_w(self.n_gen)
        self.gen_lum(self.n_gen)
        self.gen_si(self.n_gen)
        pprint(f'Finished generating {self.name} population')


if __name__ == '__main__':

    # Quick test whether everything seems to be working or not
    p = CosmicPopulation(10000)

    import matplotlib.pyplot as plt

    for arg in p.frbs.__dict__:
        print(arg)

        values = getattr(p.frbs, arg)

        if values is not None:
            plt.hist(values, bins=50)
            plt.xlabel(arg)
            plt.savefig(f'./tests/plots/{arg}.png')
            plt.clf()
