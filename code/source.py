import math


class Source:
    """Class containing individual source properties"""

    def __init__(self):

        self.flux = None
        self.dm = None
        self.dm_mw = None
        self.dm_igm = None
        self.dm_host = None
        self.width = None  # Intrinsic pulse width
        self.lum_1400 = None

        # Galactic coordinates
        self.gl = None
        self.gb = None
        self.gx = None
        self.gy = None
        self.gz = None
        self.dist = None  # Distance source to Sun [kpc]
        self.z = None

        # Detection properties
        self.snr = None
        self.w_obs = None
        self.fwhm = None
        self.s_peak = None
        self.f_obs = None
        self.detected = False

    def __str__(self):
        """Define how to print an FRB source to a console"""

        s = 'Frb source properties:'

        attributes = []
        for e in self.__dict__:
            attr = '\n\t{0:11.12}{1:.60}'.format(e, str(self.__dict__[e]))
            attributes.append(attr)

        s += ''.join(attributes)

        return s

    # Functions with which to derive properties
    def s_1400(self):
        """
        Calculate the flux density of a source at 1400 MHz

        Returns:
            s (float): Flux density of a source [mJy] at 1400 MHz
        """
        # Differs from psrpoppy in that it includes a factor 4pi
        return self.lum_1400 / (4*math.pi*self.dist**2)
