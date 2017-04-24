import random
import distributions as dis
from frb import FRB


class Source:
    """Class containing individual source properties"""

    def __init__(self):

        # "Intrinsic properties"
        self.dm = None
        self.dm_mw = None
        self.dm_igm = None
        self.dm_host = None
        self.gl = None
        self.gb = None
        self.gx = None
        self.gy = None
        self.gz = None
        self.dist = None
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

    def __str__(self):
        """Define how to print an FRB source to a console"""

        s = 'Frb source properties:'

        attributes = []
        for e in self.__dict__:
            attr = '\n\t{0:12.11}{1:.60}'.format(e, str(self.__dict__[e]))
            attributes.append(attr)

        s += ''.join(attributes)

        return s

    def add(self, frb):
        """Add an FRB to the source"""
        self.frbs.append(frb)
        self.n_frbs += 1

    def create_frb(self, pop):
        """
        Create an frb to add to source

        Args:
            pop (Population): Population parameters
        """

        # Initialise an FRB
        frb = FRB()

        # Give a redshifted random intrinsic pulse width [ms]
        frb.w_int = dis.redshift_w(z=self.z, w_min=pop.w_min, w_max=pop.w_max)

        # Add bolometric luminosity [erg/s]
        frb.lum_bol = dis.powerlaw(pop.lum_min, pop.lum_max, pop.lum_pow)

        # Add spectral index
        frb.si = random.gauss(pop.si_mean, pop.si_sigma)

        # Add FRB to source
        self.add(frb)
