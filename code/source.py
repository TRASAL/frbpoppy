import math

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


    def __str__(self):
        """Define how to print an FRB source to a console"""

        s = 'Frb source properties:'

        attributes = []
        for e in self.__dict__:
            attr = '\n\t{0:11.12}{1:.60}'.format(e, str(self.__dict__[e]))
            attributes.append(attr)

        s += ''.join(attributes)

        return s


    # From here on Galactic operations (doesn't that sound cool?!),
    # as in converting coordinates, calculating DM etc.


    def lb_to_xyz(self, dist):
        """Convert galactic coordinates to galactic XYZ"""

        rsun = 8.5  # kpc

        l = math.radians(self.gl)
        b = math.radians(self.gb)

        self.gx = dist * math.cos(b) * math.sin(l)
        self.gy = rsun - dist * math.cos(b) * math.cos(l)
        self.gz = dist * math.sin(b)
