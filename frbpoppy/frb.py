class FRB:
    """Class containing FRB properties"""

    def __init__(self):

        self.lum_bol = None
        self.si = None
        self.time = None
        self.w_arr = None
        self.w_int = None

        # Detection properties
        self.detected = False
        self.fluence = None
        self.s_peak = None
        self.snr = None
        self.w_eff = None

    def __str__(self):
        """Define how to print an FRB to a console"""

        s = 'FRB properties:'

        attributes = []
        for e in self.__dict__:
            attr = '\n\t{0:12.11}{1:.60}'.format(e, str(self.__dict__[e]))
            attributes.append(attr)

        s += ''.join(attributes)

        return s
