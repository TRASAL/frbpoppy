"""Classes to hold rate counters"""

class Rates:
    """
    Class to hold rate counters.

    Args:
        det (int, optional): Number detected
        faint (int, optional): Number to faint to detect
        jy (int, optional): Number arrived within survey > 1 Jy
        out (int, optional): Number outside survey, space or timewise
        sky (int, optional): Number per sky above > 1 Jy
        vol (int, optional): Number per Gpc^3

    """

    def __init__(self):
        """Initializing."""
        # Rates
        self.det = 0
        self.faint = 0
        self.jy = 0
        self.out = 0
        self.sky = 0
        self.vol = 0

    def tot(self):
        """Calculate the total number of rates."""
        return self.det + self.out + self.faint


class NumberOf:
    """Quick and dirty hack to store rates."""

    def __init__(self):
        """Initialise."""
        self.frbs = Rates()
        self.srcs = Rates()
