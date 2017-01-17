class Source:
    """Class containing individual source properties"""

    def __init__(self):

        self.flux = None
        self.dm = None
        self.rm = None

        # Galactic coordinates
        self.gl = None
        self.gb = None
        self.gx = None
        self.gy = None
        self.gz = None

        # Detection properties
        self.s2n = None
        self.w_obs = None
        self.fwhm = None
        self.s_peak = None
        self.f_obs = None

    def __str__(self):
        """Define how to print an FRB source to a console"""

        s = 'Frb source properties:'

        attributes = []
        for e in self.__dict__:
            attr = '\n\t{0:11.12}{1:.60}'.format(e, str(self.__dict__[e]))
            attributes.append(attr)

        s += ''.join(attributes)

        return s
