"""Class to hold FRB source properties."""
import numpy as np


class FRBs:
    """Class containing FRB properties"""

    def __init__(self):
        """Initializing."""

        # Location properties
        self.ra = None
        self.dec = None
        self.dist_co = None  # Comoving distance [Gpc]
        self.gb = None
        self.gl = None
        self.gx = None
        self.gy = None
        self.gz = None
        self.z = None

        # Dispersion measure properties
        self.dm = None
        self.dm_host = None
        self.dm_igm = None
        self.dm_mw = None

        # Intrinsic properties
        self.lum_bol = None
        self.si = None
        self.w_arr = None
        self.w_int = None

        # Detection properties
        self.fluence = None
        self.s_peak = None
        self.snr = None
        self.t_dm = 0
        self.t_dm_err = 0
        self.t_scat = 0
        self.T_sky = 0
        self.T_tot = 0
        self.w_eff = None
