"""Group together various number density descriptors."""
import random
import numpy as np
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
        m = self.dt(z=np.array([z_max]))
        self.dist_co_max = m[1]
        self.vol_co_max = m[2]
        self.cdf_sfr_max = m[-2]
        self.cdf_smd_max = m[-1]
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
        d = self.dt(vol_co=vol_co)
        z = d[0]
        dist_co = d[1]
        return z, dist_co

    def from_vol_co(self, n_gen=1):
        """Use constant number density of sources per comoving volume.

        Can be influenced by changing alpha.
        """
        vol_co = self.vol_co_max*np.random.random(n_gen)  # [Gpc]
        d = self.dt(vol_co=vol_co)
        z = d[0]
        dist_co = d[1]
        return z, dist_co

    def from_sfr(self, n_gen=1):
        """Get sources to follow star forming rate.

        Return a random redshift for sources following the Star Formation Rate.

        Follows Madau & Dickinson (2014), eq. 15. For more info see
        https://arxiv.org/pdf/1403.0007.pdf
        """
        sampling = np.random.uniform(0., self.cdf_sfr_max, size=n_gen)
        d = self.dt(cdf_sfr=sampling)
        z = d[0]
        dist_co = d[1]
        return z, dist_co

    def from_smd(self, n_gen=1):
        """
        Return a random redshift for sources following Stellar Mass Density.

        Follows Madau & Dickinson (2014), eq. 2 & 15. For more info see
        https://arxiv.org/pdf/1403.0007.pdf
        """
        sampling = np.random.uniform(0., self.cdf_smd_max, size=n_gen)
        d = self.dt(cdf_smd=sampling)
        z = d[0]
        dist_co = d[1]
        return z, dist_co
