"""Class to hold FRB source properties."""
import numpy as np
import pandas as pd
from frbpoppy.log import pprint

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

        # Repeat properties
        self.time = None

        # Detection properties
        self.fluence = None
        self.offset = None
        self.s_peak = None
        self.snr = None
        self.t_dm = 0
        self.t_scat = 0
        self.T_sky = 0
        self.T_sys = 0
        self.w_eff = None

        # Software properties
        self.index = None

    def apply(self, mask):
        """Apply a Numpy array to all parameters.

        Args:
            mask (array): Masking array to apply to all frb parameters.

        """
        for attr in self.__dict__.keys():
            parm = getattr(self, attr)
            if type(parm) is np.ndarray:
                # 1D mask on 1D array
                if mask.ndim == 1:
                    setattr(self, attr, parm[mask])
                else:
                    # 2D mask on 1D array
                    if parm.ndim == 1:
                        setattr(self, attr, parm[mask.any(axis=1)])
                    # 2D mask on 2D array
                    else:
                        parm_nan = np.where(mask, parm, np.nan)
                        row_mask = ~np.isnan(parm_nan).all(axis=1)
                        parm_nan = parm_nan[row_mask]
                        setattr(self, attr, parm_nan)

    def to_df(self):
        """Convert properties over to a Pandas DataFrame."""
        # Find all source properties
        df = pd.DataFrame()
        for attr in self.__dict__.keys():
            parm = getattr(self, attr)
            if type(parm) is np.ndarray:
                # 2D arrays to 1D
                if parm.ndim > 1:
                    df[attr] = parm[~np.isnan(self.time)]
                # 1D array to match length 2D->1D arrays
                elif type(self.time) is np.ndarray:
                    concate = np.array([parm, ]*self.time.shape[1]).transpose()
                    df[attr] = concate[~np.isnan(self.time)]
                else:
                    df[attr] = parm

        return df
