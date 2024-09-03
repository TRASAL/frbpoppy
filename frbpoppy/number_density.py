"""Group together various number density descriptors."""
import numpy as np
#import numexpr as ne
import frbpoppy.precalc as pc
import frbpoppy.galacticops as go

#global rng
#rng = np.random.default_rng()

class NumberDensity:
    """Class for cosmic number density."""

    def __init__(self,
                 model='vol_co',
                 z_max=2.0,
                 H_0=67.74,
                 W_m=0.3089,
                 W_v=0.6911,
                 alpha=-1.5,
                 delay_time=0):
        """Draw from particular number density distributions.

        Args:
            model (str): Which number density model to follow.
            z_max (float): Maximum redshift.
            H_0 (float): Hubble constant.
            W_m (float): Density parameter Ω_m.
            W_v (float): Cosmological constant Ω_Λ.
            alpha (float, optional): Desired log N log S slope for a perfect,
                non-cosmological population.
        """
        self.z_max = z_max
        self.H_0 = H_0
        self.W_m = W_m
        self.W_v = W_v
        
        # Convert various maximum distance values
        self.dt = pc.DistanceTable(H_0=H_0, W_m=W_m, W_v=W_v).lookup
        m = self.dt(z=np.array([z_max]))
        self.dist_co_max = m[1]
        self.vol_co_max = m[2]
        self.cdf_sfr_max = m[4]
        self.cdf_smd_max = m[5]
        self.cdf_dsfr0d1_max = m[6]
        self.cdf_dsfr0d5_max = m[7]
        self.cdf_dsfr1_max = m[8]
        self.dist_co_max = self.dist_co_max[0]
        self.vol_co_max = self.vol_co_max[0]
        self.alpha = alpha
        self.delay_time = delay_time
        
        # Determine from which type of distribution to draw
        if model == 'vol_co':
            self.draw = self.from_vol_co
        elif model == 'sfr':
            self.draw = self.from_sfr
        elif model == 'smd':
            self.draw = self.from_smd
        elif model == 'delayed_sfr':
            self.draw = self.from_delayed_sfr

        # Allow for the steepness of log N log S to be adapted
        if alpha != -1.5:
            self.draw = self.sloped_dist
            self.power = -self.alpha/1.5
            self.maxi = self.vol_co_max**self.power
            #self.maxi = ne.evaluate("vol_co_max**power", global_dict=vars(self))

    def sloped_dist(self, n_gen=1):
        """Draw from a sloped distribution to create a logNlogS slope."""
        
        vol_co = (self.maxi*np.random.random(n_gen))**(1/self.power)
        #rand = rng.random(n_gen, dtype=np.float32)
        #vol_co = (self.maxi * rand)**(1/self.power)  # [Gpc]
        #vol_co = ne.evaluate("(maxi*rand)**(1/power)", global_dict=vars(self))  # [Gpc]
        d = self.dt(vol_co=vol_co)
        z = d[0]
        dist_co = d[1]
        return z.astype(np.float32), dist_co.astype(np.float32)

    def from_vol_co(self, n_gen=1):
        """Use constant number density of sources per comoving volume.

        Can be influenced by changing alpha.
        """
        vol_co = self.vol_co_max*np.random.random(n_gen)
        #rand = rng.random(n_gen, dtype=np.float32)
        #vol_co = self.vol_co_max * rand  # [Gpc]
        #vol_co = ne.evaluate("vol_co_max*rand", global_dict=vars(self))  # [Gpc]
        d = self.dt(vol_co=vol_co)
        z = d[0]
        dist_co = d[1]
        return z.astype(np.float32), dist_co.astype(np.float32)

    def from_sfr(self, n_gen=1):
        """Get sources to follow star forming rate.

        Return a random redshift for sources following the Star Formation Rate.

        Follows Madau & Dickinson (2014), eq. 15. For more info see
        https://arxiv.org/pdf/1403.0007.pdf
        """
        sampling = np.random.uniform(0., self.cdf_sfr_max, size=n_gen)
        #rand = rng.random(n_gen, dtype=np.float32)
        #sampling = self.cdf_sfr_max * rand  # [Gpc]
        #sampling = ne.evaluate("cdf_sfr_max*rand", global_dict=vars(self))  # [Gpc]
        d = self.dt(cdf_sfr=sampling)
        z = d[0]
        dist_co = d[1]
        return z.astype(np.float32), dist_co.astype(np.float32)

    def from_delayed_sfr(self, n_gen=1):
        """Get sources to follow delayed star forming rate.

        Return a random redshift for sources following the delayed Star Formation Rate.

        Follows Madau & Dickinson (2014), eq. 15. For more info see
        https://arxiv.org/pdf/1403.0007.pdf
        """
        
        #rand = rng.random(n_gen, dtype=np.float32)
        
        d_list = [0.1, 0.5, 1]
        delay_time = min(d_list, key=lambda x:abs(x-self.delay_time))
        if delay_time != self.delay_time:
            print('Generate cosmic population with delay time', delay_time, 'Gyr instead')
        
        if delay_time == 0.1:
            sampling = np.random.uniform(0., self.cdf_dsfr0d1_max, size=n_gen)
            #sampling = self.cdf_dsfr0d1_max * rand
            #sampling = ne.evaluate("cdf_dsfr0d1_max*rand", global_dict=vars(self))
            d = self.dt(cdf_dsfr0d1=sampling)
        elif delay_time == 0.5:
            sampling = np.random.uniform(0., self.cdf_dsfr0d5_max, size=n_gen)
            #sampling = self.cdf_dsfr0d5_max * rand
            #sampling = ne.evaluate("cdf_dsfr0d5_max*rand", global_dict=vars(self))
            d = self.dt(cdf_dsfr0d5=sampling)
        elif delay_time == 1:
            sampling = np.random.uniform(0., self.cdf_dsfr1_max, size=n_gen)
            #sampling = self.cdf_dsfr1_max * rand
            #sampling = ne.evaluate("cdf_dsfr1_max*rand", global_dict=vars(self))
            d = self.dt(cdf_dsfr1=sampling)
        
        z = d[0]
        dist_co = d[1]
        return z.astype(np.float32), dist_co.astype(np.float32)
    
    def from_smd(self, n_gen=1):
        """
        Return a random redshift for sources following Stellar Mass Density.

        Follows Madau & Dickinson (2014), eq. 2 & 15. For more info see
        https://arxiv.org/pdf/1403.0007.pdf
        """
        sampling = np.random.uniform(0., self.cdf_smd_max, size=n_gen)
        #rand = rng.random(n_gen, dtype=np.float32)
        #sampling = self.cdf_smd_max * rand
        #sampling = ne.evaluate("cdf_smd_max*rand", global_dict=vars(self))
        d = self.dt(cdf_smd=sampling)
        z = d[0]
        dist_co = d[1]
        return z.astype(np.float32), dist_co.astype(np.float32)
