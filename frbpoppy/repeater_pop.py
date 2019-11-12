"""Allow for repeating FRB sources."""
import numpy as np
from frbpoppy.cosmic_pop import CosmicPopulation
from frbpoppy.log import pprint
from scipy.special import gamma


class RepeaterPopulation(CosmicPopulation):
    """Allow repeating FRBs to be modeled."""

    def __init__(self,
                 n_gen=1,
                 days=1,
                 generate=False,
                 lum_rep_model='independent',
                 lum_rep_sigma=1e3,
                 si_rep_model='independent',
                 si_rep_sigma=0.1,
                 times_rep_model='even',
                 w_rep_model='independent',
                 w_rep_sigma=0.05,
                 **kw):
        """Allow repeating FRB sources to be modeled.

        Args:
            lum_rep_model (str): Depenancy of bursts on prior bursts. Choice
                of 'same', 'independent', 'gaussian' or 'lognormal'.
            times_rep_model (str): 'even'-ly spaced intervals, or 'clustered'
            **kw (all): All of the arguments available to CosmicPopulation.
        """
        kw['generate'] = False
        super(RepeaterPopulation, self).__init__(n_gen, **kw)
        self.name = 'repeater'
        self.days = days
        self.n_gen = int(n_gen)

        # Time parameters
        self.times_rep_model = times_rep_model
        self.shape = n_gen

        # Luminosity parameters
        self.lum_rep_model = lum_rep_model
        self.lum_rep_sigma = lum_rep_sigma

        # Pulse width parameters
        self.w_rep_model = w_rep_model
        self.w_rep_sigma = w_rep_sigma

        # Spectral index parameters
        self.si_rep_model = si_rep_model
        self.si_rep_sigma = si_rep_sigma

        if generate:
            self.generate()

    def gen_clustered_times(self, r=5.7, k=0.34):
        """Generate burst times following Oppermann & Pen (2017)."""
        pprint('Adding burst times')
        # Determine the maximum possible number of bursts per source to include
        log_size = 1
        m = int(10**log_size)
        dims = (self.n_gen, m)

        lam = 1/(r*gamma(1 + 1/k))
        time = lam*np.random.weibull(k, dims).astype(np.float32)
        time = np.cumsum(time, axis=1)  # This is in fraction of days

        # The first burst time is actually the time since the previous one
        # You want to be at in a random time in between those
        time_offset = np.random.uniform(0, time[:, 0])
        time -= time_offset[:, np.newaxis]

        # This is interesting. The Weibull distribution was fit to an
        # observed distribution, but here I'm setting the intrinsic one
        if isinstance(self.frbs.z, np.ndarray):
            time = time*(1+self.frbs.z)[:, np.newaxis]
        else:
            time = time*(1+self.frbs.z)

        # Mask any frbs over the maximum time (Earth perspective)
        time[(time > self.days)] = np.nan

        # Iteratively add extra bursts until over limit
        mask = ~np.isnan(time[:, -1])  # Where more bursts are needed
        sum_mask = np.count_nonzero(mask)
        while sum_mask != 0:  # Add additional columns
            m = int(10**log_size)
            ms = np.full((self.n_gen, m), np.nan)
            new = lam*np.random.weibull(k, (sum_mask, m)).astype(np.float32)
            new = np.cumsum(new, axis=1)

            # Add redshift correction
            z = self.frbs.z
            if isinstance(self.frbs.z, np.ndarray):
                z = self.frbs.z[mask][:, np.newaxis]
            new *= (1+z)

            new += time[:, -1][mask][:, np.newaxis]  # Ensure cumulative
            new[(new > self.days)] = np.nan  # Apply filter
            ms[mask] = new  # Set up additional columns
            time = np.hstack((time, ms))  # Add to original array

            # Set up for next loop
            if sum_mask == self.n_gen:
                log_size += 0.5
            mask = ~np.isnan(time[:, -1])
            sum_mask = np.count_nonzero(mask)

        # Remove any columns filled with NaNs
        time = time[:, ~np.all(np.isnan(time), axis=0)]

        self.frbs.time = time

    def gen_regular_times(self):
        """Generate a series of regular spaced times."""
        n = 6  # n bursts per day
        u = np.random.uniform
        time = u(0, self.days, (self.n_gen, n)).astype(np.float32)
        time = np.sort(time)
        self.frbs.time = time*(1+self.frbs.z)[:, np.newaxis]

    def gen_rep_times(self):
        """Generate repetition times."""
        if self.times_rep_model == 'clustered':
            self.gen_clustered_times()
        else:
            self.gen_regular_times()

        # Set size for all other parameters
        self.shape = self.frbs.time.shape

    def gen_rep_w(self):
        """Add extra pulse width values."""
        if self.w_rep_model == 'same':
            self.gen_w(self.n_gen)

        elif self.w_rep_model == 'independent':
            self.gen_w(self.shape)

        elif self.w_rep_model == 'gaussian':
            self.gen_w(self.n_gen)
            mu = self.frbs.w_int
            sigma = self.w_rep_sigma
            shape = (self.shape[-1], self.shape[0])
            self.frbs.w_int = np.random.normal(mu, sigma, shape).T
            self.frbs.w_int.astype(np.float32)
            # Recalculate the pulse width upon arrival to Earth
            self.frbs.w_arr = self.frbs.w_int*(1+self.frbs.z[:, None])

    def gen_rep_si(self):
        """Generate spectral indices for repeat bursts."""
        if self.si_rep_model == 'same':
            self.gen_si(self.n_gen)

        elif self.si_rep_model == 'independent':
            self.gen_si(self.shape)

        elif self.si_rep_model == 'gaussian':
            self.gen_si(self.n_gen)
            mu = self.frbs.si
            sigma = self.si_rep_sigma
            shape = (self.shape[-1], self.shape[0])
            self.frbs.si = np.random.normal(mu, sigma, shape).T
            self.frbs.si.astype(np.float64)

    def gen_rep_lum(self):
        """Generate luminosities for repeat bursts."""
        if self.lum_rep_model == 'same':
            self.gen_lum(self.n_gen)

        elif self.lum_rep_model == 'independent':
            self.gen_lum(self.shape)

        elif self.lum_rep_model in ('gaussian', 'lognormal'):
            if self.lum_rep_model == 'gaussian':
                r = np.random.normal
            else:
                r = np.random.lognormal

            self.gen_lum(self.n_gen)
            mu = self.frbs.lum_bol
            sigma = self.lum_rep_sigma
            shape = (self.shape[-1], self.shape[0])
            self.frbs.lum_bol = r(mu, sigma, shape).T
            self.frbs.lum_bol.astype(np.float64)

    def generate(self):
        """Generate a Repeater Population."""
        pprint(f'Generating {self.name} population')
        self.gen_dist()
        self.gen_direction()
        self.gen_gal_coords()
        self.gen_dm()
        self.gen_rep_times()
        self.gen_rep_lum()
        self.gen_rep_w()
        self.gen_rep_si()
        pprint(f'Finished generating {self.name} population')

    @classmethod
    def simple(cls, n, generate=False):
        """Set up a simple, local population."""
        pop = cls(n,
                  days=1,
                  name='simple',
                  H_0=67.74,
                  W_m=0.3089,
                  W_v=0.6911,
                  dm_host_model='gaussian',
                  dm_host_mu=0.,
                  dm_host_sigma=0.,
                  dm_igm_index=0.,
                  dm_igm_sigma=None,
                  dm_mw_model='zero',
                  emission_range=[10e6, 10e9],
                  lum_range=[1e38, 1e38],
                  lum_index=0.,
                  n_model='vol_co',
                  alpha=-1.5,
                  w_model='uniform',
                  w_range=[10, 10],
                  w_mu=0.1,
                  w_sigma=1.,
                  si_mu=0.,
                  si_sigma=0.,
                  z_max=0.01,
                  lum_rep_model='independent',
                  lum_rep_sigma=1e3,
                  si_rep_model='same',
                  si_rep_sigma=0.1,
                  times_rep_model='even',
                  w_rep_model='independent',
                  generate=generate)
        return pop


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    pop = RepeaterPopulation.simple(1e5)
    pop.frbs.z = 0
    pop.gen_clustered_times(r=10, k=0.3)
    n_bursts = (~np.isnan(pop.frbs.time)).sum(1)
    plt.hist(n_bursts, bins=max(n_bursts))
    plt.yscale('log')
    plt.show()
