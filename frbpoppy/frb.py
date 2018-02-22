class FRB:
    """Class containing FRB properties"""

    def __init__(self):

        self.w_int = None
        self.lum_bol = None
        self.si = None
        self.time = None

        # Detection properties
        self.snr = None
        self.detected = False
        self.w_eff = None
        self.s_peak = None
        self.fluence = None

    def __str__(self):
        """Define how to print an FRB to a console"""

        s = 'FRB properties:'

        attributes = []
        for e in self.__dict__:
            attr = '\n\t{0:12.11}{1:.60}'.format(e, str(self.__dict__[e]))
            attributes.append(attr)

        s += ''.join(attributes)

        return s