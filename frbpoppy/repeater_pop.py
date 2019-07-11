"""Allow for repeating FRB sources."""
import numpy as np
from frbpoppy.cosmic_pop import CosmicPopulation
from frbpoppy.log import pprint


class RepeaterPopulation(CosmicPopulation):
    """Allow repeating FRBs to be modeled."""

    def __init__(self,
                 n_gen=1,
                 generate=False,
                 **kw):
        """Allow repeating FRB sources to be modeled.

        Args:
            frac (float): Fraction of population to be repeaters.
            **kw (all): All of the arguments available to CosmicPopulation.
        """
        kw['generate'] = False
        super(RepeaterPopulation, self).__init__(n_gen, **kw)
        self.name = 'repeater'
        self.n_gen = n_gen

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
            # Ensure units are per second rather than per day
            time = 86400*r*np.random.weibull(k, (n, m)).astype(np.float32)
            time = np.cumsum(time, axis=1)
            return time

        def accept(time, max_time=43200):
            """Set a limit on the maximum total bursting time.

            Args:
                times (type): Burst times
                max_time (type): Maximum of 12 hours (12*60*60 seconds)

            Returns:
                array: Mask with True if at least one burst under max time

            """
            cond = time < max_time
            time[~cond] = np.nan
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

        self.frbs.bursts = ~np.isnan(time)
        self.frbs.time = time

    def gen_regular_times(self):
        """Generate a series of regular spaced times."""
        # 24 bursts per day
        day = np.linspace(0, 1, 24).astype(np.float32)
        time = np.tile(day, (self.n_gen, 1))
        self.frbs.time = time*(1+self.frbs.z)[:, np.newaxis]

    def gen_dm(self):
        """Adapt initial DM values for subsequent bursts."""
        pass

    def gen_w_eff(self):
        """Add extra pulse width values."""
        pass

    def generate(self):
        """Gen Cosmic Population then transform into RepeaterPopulation."""
        super(RepeaterPopulation, self).generate()
        # self.gen_clustered_times()
        self.gen_regular_times()


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    n = 2000
    pop = RepeaterPopulation(n, generate=True)
    # n_bursts = (~np.isnan(pop.frbs.time)).sum(1)
    # plt.hist(n_bursts, bins=[i for i in range(11)])
    # plt.yscale('log')
    # plt.show()
    # print(pop.n_miss, pop.n_miss/n)
