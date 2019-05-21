"""Class to hold FRB source properties."""
import numpy as np
import pandas as pd

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
        self.offset = None
        self.s_peak = None
        self.snr = None
        self.t_dm = 0
        self.t_scat = 0
        self.T_sky = 0
        self.T_tot = 0
        self.w_eff = None

    def apply(self, mask):
        """Apply a Numpy array to all parameters.

        Args:
            mask (array): Masking array to apply to all frb parameters.

        """
        for attr in self.__dict__.keys():
            parm = getattr(self, attr)
            if type(parm) is np.ndarray:
                setattr(self, attr, parm[mask])

    def to_df(self):
        """Convert properties over to a Pandas DataFrame."""
        # Find all source properties
        df = pd.DataFrame()
        for attr in self.__dict__.keys():
            parm = getattr(self, attr)
            if type(parm) is np.ndarray:
                df[attr] = parm

        return df
