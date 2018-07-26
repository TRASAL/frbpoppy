"""Class containing properties relevant to an FRB source."""
import random
import frbpoppy.distributions as dis
from frbpoppy.frb import FRB


class Source:
    """Class containing individual source properties."""

    def __init__(self):
        """Initializing."""
        # "Intrinsic properties"
        self.dec = None
        self.dist = None
        self.dist_co = None  # Comoving distance [Gpc]
        self.dm = None
        self.dm_host = None
        self.dm_igm = None
        self.dm_mw = None
        self.gb = None
        self.gl = None
        self.gx = None
        self.gy = None
        self.gz = None
        self.ra = None
        self.z = None

        # Collect all FRB bursts
        self.frbs = []
        self.n_frbs = 0

        # Observing properties
        self.t_dm = 0
        self.t_dm_err = 0
        self.t_scat = 0
        self.T_sky = 0
        self.T_tot = 0
        self.detected = False

        # Special property
        self.name = None

    def __str__(self):
        """Define how to print an FRB source to a console."""
        s = 'Frb source properties:'

        attributes = []
        for e in self.__dict__:
            attr = '\n\t{0:12.11}{1:.60}'.format(e, str(self.__dict__[e]))
            attributes.append(attr)

        s += ''.join(attributes)

        return s

    def add(self, frb):
        """Add an FRB to the source."""
        self.frbs.append(frb)
        self.n_frbs += 1

    def create_frb(self, pop, time=None):
        """
        Create an frb to add to source.

        Args:
            pop (Population): Population parameters
            time (float): Time of burst [s]
        """
        # Initialise an FRB
        frb = FRB()

        # Get a random intrinsic pulse width [ms]
        if pop.w_model == 'lognormal':
            frb.w_int = random.lognormvariate(pop.w_mu, pop.w_sigma)

        if pop.w_model == 'uniform':
            frb.w_int = random.uniform(pop.w_min, pop.w_max)

        # Calculate the pulse width upon arrival to Earth
        frb.w_arr = frb.w_int*(1+self.z)

        # Add bolometric luminosity [erg/s]
        frb.lum_bol = dis.powerlaw(pop.lum_min, pop.lum_max, pop.lum_pow)

        # Add spectral index
        frb.si = random.gauss(pop.si_mu, pop.si_sigma)

        # Add time
        frb.time = time

        # Add FRB to source
        self.add(frb)
