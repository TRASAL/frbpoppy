"""Group together various number density descriptors."""
import random
import numpy as np
import frbpoppy.distributions as dis
import frbpoppy.precalc as pc


class NumberDensity:
    """Class for cosmic number density."""

    def __init__(self,
                 model='vol_co',
                 z_max=6.0,
                 H_0=67.74,
                 W_m=0.3089,
                 W_v=0.6911,
                 alpha=-1.5):
        """Draw from particular number density distributions.

        Args:
            model (str): Which number density to follow.
            z_max (float): Maximum redshift to which to draw.
            vol_co_max (float, optional): Maximum comoving redshift [Gpc^3].
            dist_co_max (float, optional): Maximum comoving distance [Gpc].
            alpha (float, optional): Desired log N log S slope for a perfect,
                non-cosmological population.
        """
        self.z_max = z_max

        # Convert various maximum distance values
        self.dt = pc.DistanceTable(H_0=H_0, W_m=W_m, W_v=W_v).lookup
        _, self.dist_co_max, self.vol_co_max = self.dt(z=np.array([z_max]))
        self.dist_co_max = self.dist_co_max[0]
        self.vol_co_max = self.vol_co_max[0]
        self.alpha = alpha

        # Determine from which type of distribution to draw
        if model == 'vol_co':
            self.draw = self.from_vol_co
        elif model == 'sfr':
            self.draw = self.from_sfr
        elif model == 'smd':
            self.draw = self.from_smd

        # Allow for the steepness of log N log S to be adapted
        if alpha != -1.5:
            self.draw = self.sloped_dist
            self.power = -self.alpha/1.5
            self.maxi = self.vol_co_max**self.power

    def sloped_dist(self, n_gen=1):
        """Draw from a sloped distribution to create a logNlogS slope."""
        vol_co = (self.maxi*np.random.random(n_gen))**(1/self.power)  # [Gpc]
        z, dist_co, _ = self.dt(vol_co=vol_co)
        return z, dist_co

    def from_vol_co(self, n_gen=1):
        """Use constant number density of sources per comoving volume.

        Can be influenced by changing alpha.
        """
        vol_co = self.vol_co_max*np.random.random(n_gen)  # [Gpc]
        z, dist_co, _ = self.dt(vol_co=vol_co)
        return z, dist_co

    def from_sfr(self, n_gen=1):
        """Get sources to follow star forming rate."""
        z = dis.z_from_sfr(z_max=self.z_max, n_gen=n_gen)
        _, dist_co, _ = self.dt(z=z)
        return z, dist_co

    def from_smd(self, n_gen=1):
        """Get sources to follow stellar mass density."""
        z = dis.z_from_csmd(z_max=self.z_max, n_gen=n_gen)
        _, dist_co, _ = self.dt(z=z)
        return z, dist_co
