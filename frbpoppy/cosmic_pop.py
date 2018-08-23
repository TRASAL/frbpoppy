"""Class to generate a cosmic population of FRBs."""
import math
import random

from frbpoppy.log import pprint
from frbpoppy.population import Population
from frbpoppy.source import Source
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
                 dm_host=100,
                 dm_igm_index=1200,
                 dm_mw_model='ne2001',
                 emission_range=[10e6, 10e9],
                 lum_range=[1e40, 1e50],
                 lum_index=0,
                 n_model='sfr',
                 pulse_model='lognormal',
                 pulse_range=[0.1, 10],
                 pulse_mu=1.6,
                 pulse_sigma=1.,
                 repeat=0.0,
                 si_mu=-1.4,
                 si_sigma=1.,
                 z_max=2.5,
                 test=False):
        """Generate a popuation of FRBs.

        Args:
            n_gen (int): Number of FRB sources/sky/time to generate.
            days (float): Number of days over which FRBs are generated.
            name (str): Population name.
            H_0 (float): Hubble constant.
            W_m (float): Density parameter Ω_m.
            W_v (float): Cosmological constant Ω_Λ.
            dm_host (float): Dispersion measure host [pc/cm^3].
            dm_igm_index (float): Dispersion measure slope for IGM [pc/cm^3].
            dm_mw_model (str): Dispersion measure model for the Milky Way.
                Options are 'ne2001' or 'zero'.
            emission_range (list): The frequency range [Hz] between which FRB
                sources should emit the given bolometric luminosity.
            lum_range (list): Bolometric luminosity (distance) range [erg/s].
            lum_index (float): Power law index.
            n_model (str): Number density model. Either 'constant' or 'sfr'.
            pulse_model (str): Pulse width model, 'lognormal' or 'uniform'.
            pulse_range (list): Pulse width range [ms].
            pulse_mu (float): Mean pulse width [ms].
            pulse_sigma (float): Deviation pulse width [ms].
            repeat (float): Repeater fraction of population.
            si_mu (float): Mean spectral index.
            si_sigma (float): Standard deviation spectral index.
            z_max (float): Maximum redshift.
            test (bool): Testing flag.

        Returns:
            Population: Population of FRBs.

        """
        # Set up population
        Population.__init__(self)
        self.dm_host = dm_host
        self.dm_igm_index = dm_igm_index
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
        self.repeat = repeat
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
        r = go.Redshift(self.z_max, H_0=self.H_0, W_m=self.W_m, W_v=self.W_v)
        self.dist_co_max = r.dist_co()
        self.vol_co_max = r.vol_co()

        # Ensure precalculations are done if necessary
        pc.dist_table(z=1.0,
                      dist_co=True,
                      H_0=self.H_0,
                      W_m=self.W_m,
                      W_v=self.W_v)

        # Let user know what's happening
        pprint(f'Generating {self.name} population')

        while self.n_srcs < self.n_gen:

            # Initialise
            src = Source()

            # Add random directional coordinates
            src.gl = random.random() * 360.0 - 180
            src.gb = math.degrees(math.asin(random.random()))
            if random.random() < 0.5:
                src.gb *= -1

            # Convert
            src.ra, src.dec = go.lb_to_radec(src.gl, src.gb)

            # Use constant number density of sources per comoving volume
            if self.n_model == 'constant':
                # Calculate comoving volume [Gpc]
                vol_co = self.vol_co_max*random.random()
                r = pc.dist_table(vol_co=vol_co, z=True, dist_co=True)
                src.dist_co = r['dist_co']
                src.z = r['z']

            # Get sources to follow star forming rate
            if self.n_model == 'sfr':
                src.z = dis.z_from_sfr(z_max=self.z_max)
                src.dist_co = pc.dist_table(z=src.z, dist_co=True)

            # Get the proper distance
            dist_pr = src.dist_co/(1+src.z)

            # Convert into galactic coordinates
            src.gx, src.gy, src.gz = go.lb_to_xyz(src.gl, src.gb, dist_pr)

            # Dispersion measure of the Milky Way
            if self.dm_mw_model == 'ne2001':
                src.dm_mw = pc.ne2001_table(src.gl, src.gb)

            elif self.dm_mw_model == 'zero':
                src.dm_mw = 0.

            # Dispersion measure of the intergalactic medium
            src.dm_igm = go.ioka_dm_igm(src.z, slope=self.dm_igm_index)

            # Dispersion measure of the host (Tendulkar)
            src.dm_host = self.dm_host / (1 + src.z)

            # Total dispersion measure
            src.dm = src.dm_mw + src.dm_igm + src.dm_host

            # Add initial frb
            attrs = dict(self.__dict__)
            attrs.pop('time')
            src.create_frb(**attrs)

            # If repeating add another FRB
            if random.random() < self.repeat:

                times = dis.oppermann_pen()
                for t in times:
                    src.create_frb(**attrs, time=t)

            # Add source to population
            self.add(src)

        pprint(f'Finished')
