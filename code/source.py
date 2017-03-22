import math


class Source:
    """Class containing individual source properties"""

    def __init__(self):

        self.flux = None
        self.dm = None
        self.dm_mw = None
        self.dm_igm = None
        self.dm_host = None
        self.w_int = None  # Intrinsic pulse width [ms]
        self.lum_bol = None
        self.si = None  # Spectral index

        # Galactic coordinates
        self.gl = None
        self.gb = None
        self.gx = None
        self.gy = None
        self.gz = None
        self.dist = None  # Distance source to Sun [Gpc]
        self.z = None

        # Detection properties
        self.snr = None
        self.fwhm = None
        self.detected = False
        self.w_eff = None
        self.s_peak = None
        self.fluence = None

    def __str__(self):
        """Define how to print an FRB source to a console"""

        s = 'Frb source properties:'

        attributes = []
        for e in self.__dict__:
            attr = '\n\t{0:12.11}{1:.60}'.format(e, str(self.__dict__[e]))
            attributes.append(attr)

        s += ''.join(attributes)

        return s
