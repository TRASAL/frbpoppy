"""Group together various number density descriptors."""
import random
import math
import frbpoppy.distributions as dis
import frbpoppy.precalc as pc


class NumberDensity:
    """Class for cosmic number density."""

    def __init__(self,
                 model='vol_co',
                 z_max=8.0,
                 vol_co_max=3201.21717511309,
                 dist_co_max=9.14295802075115,
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
        self.vol_co_max = vol_co_max
        self.dist_co_max = dist_co_max
        self.alpha = alpha

        if model == 'vol_co':
            self.draw = self.from_vol_co
        elif model == 'sfr':
            self.draw = self.from_sfr
        elif model == 'smd':
            self.draw = self.from_smd

        if alpha == -1.5:
            self.random_num = random.random
        else:
            # See sloped_dist
            # self.a = 2*math.log(-self.alpha/3)/math.log(2) + 2  # Slope of PDF
            # self.b = 0.5-(self.a/2)  # Offset of PDF
            self.random_num = self.sloped_dist
            self.m = -2*self.alpha - 3

    def sloped_dist(self):
        """Draw from a sloped distribution to create a certain logNlogS slope.

        For PDF=a*x+b, and CDF=c*((a/2)*x**2+b*x) with c=2 so that P(.5)=.5
        """
        y = random.random()  # Uniform distribution
        # return (-self.b + math.sqrt(self.a*y+self.b**2))/self.a
        return dis.powerlaw(low=1e-99, high=1, power=self.m)

    def from_vol_co(self):
        """Use constant number density of sources per comoving volume.

        Can be influenced by changing alpha.
        """
        vol_co = self.vol_co_max*self.random_num()  # [Gpc]
        r = pc.dist_table(vol_co=vol_co, z=True, dist_co=True)
        dist_co = r['dist_co']
        z = r['z']
        return dist_co, z

    def from_sfr(self):
        """Get sources to follow star forming rate."""
        z = dis.z_from_sfr(z_max=self.z_max)
        dist_co = pc.dist_table(z=z, dist_co=True)
        return dist_co, z

    def from_smd(self):
        """Get sources to follow stellar mass density."""
        z = dis.z_from_csmd(z_max=self.z_max)
        dist_co = pc.dist_table(z=z, dist_co=True)
        return dist_co, z
