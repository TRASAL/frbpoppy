# -*- coding: future_fstrings -*-
"""Allow for repeating FRB sources."""
import numpy as np
from frbpoppy.cosmic_pop import CosmicPopulation
from frbpoppy.log import pprint


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
                of 'same', 'independent', 'normal' or 'lognormal'.
            times_rep_model (str): 'even'-ly spaced intervals, or 'clustered'
            **kw (all): All of the arguments available to CosmicPopulation.
        """
        kw['generate'] = False
        super(RepeaterPopulation, self).__init__(n_gen, **kw)
        self.name = 'repeater'
        self.days = days
        self.n_gen = n_gen

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

    def gen_clustered_times(self):
        """Generate burst times following Oppermann & Pen (2017).

        As most bursts will be beyond 12h, it uses a accept-reject algorithm to
        ensure all frbs are useful in the surveys phase. It uses a counter to
        keep track of the number of FRBs it has had to regenerate.
        """
        pprint('Adding burst times')
        # Define parameters for Weibull distribution
        r = 5.7
        k = 0.34

        # Determine the maximum possible number of bursts to include
        # Would be minus one if bursts over 12h were included
        m = int(round(np.log10(self.n_gen))*2)
        if m < 1:
            m = 1
        dims = (self.n_gen, m)

        def sample(n, m):
            # Ensure unit per day
            time = r*np.random.weibull(k, (n, m)).astype(np.float32)
            time = np.cumsum(time, axis=1)
            return time

        def accept(time, max_time=1.):
            """Set a limit on the maximum total bursting time.

            Args:
                times (type): Burst times
                max_time (type): Maximum time to still accept

            Returns:
                array: Mask with True if at least one burst under max time

            """
            time[time > max_time] = np.nan
            return ~np.isnan(time[:, 0])

        # Accept - reject algorithm
        self.n_miss = 0
        time = sample(*dims)
        mask = accept(time)
        reject, = np.where(~mask)
        self.n_miss = reject.size
        while reject.size > 0:
            fill = sample(reject.size, m)
            mask = accept(fill)
            time[reject[mask]] = fill[mask]
            reject = reject[~mask]

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

        elif self.w_rep_model == 'normal':
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

        elif self.si_rep_model == 'normal':
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

        elif self.lum_rep_model in ('normal', 'lognormal'):
            if self.lum_rep_model == 'normal':
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


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    n = 2000
    pop = RepeaterPopulation(n, times_rep_model='clustered', generate=True)
    n_bursts = (~np.isnan(pop.frbs.time)).sum(1)
    print(n_bursts)
    plt.hist(n_bursts, bins=10)
    plt.yscale('log')
    # plt.show()
