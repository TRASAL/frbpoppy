"""Allow for repeating FRB sources."""
import numpy as np
from frbpoppy.cosmic_pop import CosmicPopulation
from frbpoppy.log import pprint
import frbpoppy.time_dists as td
import frbpoppy.w_dists as wd
import frbpoppy.si_dists as sid
import frbpoppy.lum_dists as ld


class RepeaterPopulation(CosmicPopulation):
    """Allow repeating FRBs to be modeled."""

    def __init__(self,
                 n_srcs=1,
                 n_days=1,
                 repeaters=False,
                 generate=False,
                 **kw):
        """Allow repeating FRB sources to be modeled.

        Args:
            **kw (all): All of the arguments available to CosmicPopulation.
        """
        kw['generate'] = False
        super(RepeaterPopulation, self).__init__(n_srcs, **kw)
        self.repeaters = repeaters
        self.n_days = n_days
        self.n_srcs = int(n_srcs)
        self.shape = (self.n_srcs,)

        if self.repeaters:
            self.name = 'repeater'
            self.set_time()

        self.set_lum()
        self.set_si()
        self.set_w()

        if generate:
            self.generate()

    def set_time(self, f='regular', **kwargs):
        """Set function from which to generate time stamps.

        Usage is in the form of `pop.set_time('clustered', r=2, k=5)`, and then
        then parameters can be generated using pop.gen_time()

        Options:
            'regular'
                Args:
                    lam: Number of bursts per day
            'poisson'
                Args:
                    lam: Expected number of bursts per day
            'clustered'
                Args:
                    r: Rate parameter
                    k: Shape parameter
        Args:
            f (str/function): Either string to link to function or function
                itself to use for generating parameter.
            **kwargs (dict): Keywords to pass on to corresponding function.

        """
        if not isinstance(f, str):
            # These lambda functions look complex, but aren't.
            # They merely stop the function from running immediately
            self.time_func = lambda: f(**kwargs)
            return

        if f == 'regular':
            self.time_func = lambda: td.regular(n_srcs=self.n_srcs,
                                                n_days=self.n_days,
                                                z=self.frbs.z,
                                                **kwargs)
        elif f.startswith('poisson'):
            self.time_func = lambda: td.poisson(n_srcs=self.n_srcs,
                                                n_days=self.n_days,
                                                z=self.frbs.z,
                                                **kwargs)
        elif f == 'clustered':
            self.time_func = lambda: td.clustered(n_srcs=self.n_srcs,
                                                  n_days=self.n_days,
                                                  z=self.frbs.z,
                                                  **kwargs)
        else:
            raise IOError('set_time input not recognised')

    def gen_time(self):
        """Generate times."""
        self.frbs.time = self.time_func()
        # Set size for all other parameters
        self.shape = self.frbs.time.shape

    def set_w(self, f='uniform', per_source='same', **kwargs):
        """Set pulse widths functions [ms]."""
        if not isinstance(f, str):
            self.time_func = lambda: f(**kwargs)
            return

        # Each burst from the same source: same or different widths?
        if per_source == 'same':
            shape = self.n_srcs
        elif per_source == 'different':
            shape = self.shape
        # Or use a gaussian per source.
        elif per_source.startswith('gauss'):
            self.w_func = lambda: wd.gauss_per_source(dist=getattr(wd, f),
                                                      shape=self.shape,
                                                      z=self.frbs.z,
                                                      **kwargs)
            return

        # Distribution from which to draw pulse widths
        if f == 'uniform':
            self.w_func = lambda: wd.uniform(shape=shape, z=self.frbs.z,
                                             **kwargs)
        elif f == 'lognormal':
            self.w_func = lambda: wd.lognormal(shape=shape, z=self.frbs.z,
                                               **kwargs)
        else:
            raise IOError('set_w input not recognised')

    def gen_w(self):
        """Generate pulse widths [ms]."""
        self.frbs.w_int, self.frbs.w_arr = self.w_func()

    def set_si(self, f='gauss', per_source='same', **kwargs):
        """Set spectral index functions [ms]."""
        if not isinstance(f, str):
            self.time_func = lambda: f(**kwargs)
            return

        # Each burst from the same source: same or different siidths?
        if per_source == 'same':
            shape = self.n_srcs
        elif per_source == 'different':
            shape = self.shape
        # Or use a gaussian per source
        elif per_source.startswith('gauss'):
            dist = getattr(sid, f)
            self.si_func = lambda: sid.gauss_per_source(dist=dist,
                                                        shape=self.shape,
                                                        **kwargs)
            return

        # Distribution from sihich to drasi pulse siidths
        if f.startswith('gauss'):
            self.si_func = lambda: sid.gauss(shape=shape, **kwargs)
        else:
            raise IOError('set_si input not recognised')

    def gen_si(self):
        """Generate spectral indices."""
        self.frbs.si = self.si_func()

    def set_lum(self, f='powerlaw', per_source='same', **kwargs):
        """Set luminosity function [ergs/s]."""
        if not isinstance(f, str):
            self.time_func = lambda: f(**kwargs)
            return

        # Each burst from the same source: same or different luminosities?
        if per_source == 'same':
            shape = self.n_srcs
        elif per_source == 'different':
            shape = self.shape
        # Or use a gaussian per source
        elif per_source.startswith('gauss'):
            dist = getattr(ld, f)
            self.lum_func = lambda: ld.gauss_per_source(dist=dist,
                                                        shape=self.shape,
                                                        **kwargs)
            return

        if f == 'powerlaw':
            self.lum_func = lambda: ld.powerlaw(shape=shape, **kwargs)
        else:
            raise IOError('set_lum input not recognised')

    def gen_lum(self):
        """Generate luminosities [ergs/s]."""
        self.frbs.lum_bol = self.lum_func()

    def generate(self):
        """Generate a Repeater Population."""
        pprint(f'Generating {self.name} population')
        self.gen_dist()
        self.gen_direction()
        self.gen_gal_coords()
        self.gen_dm()
        self.gen_time()
        self.gen_lum()
        self.gen_si()
        self.gen_w()
        pprint(f'Finished generating {self.name} population')

    @classmethod
    def simple(cls, n, generate=False):
        """Set up a simple, local population."""
        # # TODO: UPDATE!
        pop = cls(n,
                  n_days=1,
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
                  emission_range=[100e6, 10e9],
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
                  lum_rep_model='same',
                  lum_rep_sigma=1e3,
                  si_rep_model='same',
                  si_rep_sigma=0.1,
                  w_rep_model='same',
                  generate=generate)
        pop.set_time('regular')
        return pop


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    pop = RepeaterPopulation.simple(1e4)
    pop.generate()
    import IPython; IPython.embed()
    # plt.hist(pop.frbs.si)
    # plt.show()
    # n_bursts = (~np.isnan(pop.frbs.time)).sum(1)
    # plt.hist(n_bursts, bins=max(n_bursts))
    # plt.yscale('log')
    # plt.show()
