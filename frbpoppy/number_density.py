"""Group together various number density descriptors."""
import random


class NumberDensity:

    def __init__(self,
                 z_max=8.0,
                 vol_co_max=3201.21717511309,
                 dist_co_max=9.14295802075115):
        """TODO"""
        self.z_max = z_max
        self.vol_co_max = vol_co_max
        self.dist_co_max = dist_co_max

    def from_vol_co(self, alpha=None):
        """Use constant number density of sources per comoving volume."""
        vol_co = self.vol_co_max*random.random()  # [Gpc]
        r = pc.dist_table(vol_co=vol_co, z=True, dist_co=True)
        dist_co = r['dist_co']
        z = r['z']
        return dist_co, z

    def from_sfr(self):
        """Get sources to follow star forming rate."""
        z = dis.z_from_sfr(z_max=self.z_max)
        dist_co = pc.dist_table(z=z, dist_co=True)

    def from_smd(self):
        """Get sources to follow stellar mass density."""
        z = dis.z_from_csmd(z_max=self.z_max)
        dist_co = pc.dist_table(z=z, dist_co=True)

    def from_z(self):
        """Get sources to follow redshift."""
        z = self.z_max*random.random()
        r = pc.dist_table(z=z, dist_co=True, vol_co=True)
        dist_co = r['dist_co']
        vol_co = r['vol_co']

    def from_dist_co(self):
        """Get sources to follow comoving distance"""
        dist_co = self.dist_co_max*random.random()
        r = pc.dist_table(dist_co=dist_co, z=True, vol_co=True)
        z = r['z']
        vol_co = r['vol_co']
