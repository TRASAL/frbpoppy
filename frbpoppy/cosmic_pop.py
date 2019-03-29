"""Class to generate a cosmic population of FRBs."""
import numpy as np
import math
import random

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
                 H_0=69.6,
                 W_m=0.286,
                 W_v=0.714,
                 dm_host_model='normal',
                 dm_host_mu=100,
                 dm_host_sigma=0,
                 dm_igm_index=1000,
                 dm_igm_sigma=None,
                 dm_mw_model='ne2001',
                 emission_range=[10e6, 10e9],
                 lum_range=[1e40, 1e45],
                 lum_index=0,
                 n_model='sfr',
                 alpha=-1.5,
                 pulse_model='lognormal',
                 pulse_range=[0.1, 10],
                 pulse_mu=1.6,
                 pulse_sigma=1.,
                 si_mu=-1.4,
                 si_sigma=1.,
                 z_max=2.5):
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
            pulse_model (str): Pulse width model, 'lognormal' or 'uniform'.
            pulse_range (list): Pulse width range [ms].
            pulse_mu (float): Mean pulse width [ms].
            pulse_sigma (float): Deviation pulse width [ms].
            si_mu (float): Mean spectral index.
            si_sigma (float): Standard deviation spectral index.
            z_max (float): Maximum redshift.

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
        self.n_gen = n_gen
        self.n_model = n_model
        self.si_mu = si_mu
        self.si_sigma = si_sigma
        self.time = days * 86400  # Convert to seconds
        self.w_model = pulse_model
        self.w_max = pulse_range[1]
        self.w_min = pulse_range[0]
        self.w_mu = pulse_mu
        self.w_sigma = pulse_sigma
        self.W_m = W_m
        self.W_v = W_v
        self.z_max = z_max

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

        # Let user know what's happening
        pprint(f'Generating {self.name} population')

        frbs = self.frbs

        # Add random directional coordinates
        frbs.gl = np.random.random(n_gen) * 360.0 - 180
        frbs.gb = np.degrees(np.arcsin(np.random.random(n_gen)))
        frbs.gb[::2] *= -1

        # Convert
        frbs.ra, frbs.dec = go.lb_to_radec(frbs.gl, frbs.gb)

        # Draw from number density
        frbs.z, frbs.dist_co = n_den(n_gen)

        # Get the proper distance
        dist_pr = frbs.dist_co/(1+frbs.z)

        # Convert into galactic coordinates
        frbs.gx, frbs.gy, frbs.gz = go.lb_to_xyz(frbs.gl, frbs.gb, dist_pr)

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
        if self.dm_host_model == 'normal':
            frbs.dm_host = np.random.normal(self.dm_host_mu,
                                            self.dm_host_sigma)
        elif self.dm_host_model == 'lognormal':
            mu = math.log(self.dm_host_mu)
            sigma = math.log(self.dm_host_sigma)
            frbs.dm_host = np.random.lognormal(mu, sigma)

        frbs.dm_host /= (1 + frbs.z)

        # Total dispersion measure
        frbs.dm = frbs.dm_mw + frbs.dm_igm + frbs.dm_host

        # Get a random intrinsic pulse width [ms]
        if self.w_model == 'lognormal':
            frbs.w_int = np.random.lognormal(self.w_mu, self.w_sigma, n_gen)

        if self.w_model == 'uniform':
            frbs.w_int = np.random.uniform(self.w_min, self.w_max, n_gen)

        # Calculate the pulse width upon arrival to Earth
        frbs.w_arr = frbs.w_int*(1+frbs.z)

        # Add bolometric luminosity [erg/s]
        frbs.lum_bol = dis.powerlaw(self.lum_min,
                                    self.lum_max,
                                    self.lum_pow,
                                    n_gen)

        # Add spectral index
        frbs.si = np.random.normal(si_mu, si_sigma, n_gen)

        pprint('Finished')


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
