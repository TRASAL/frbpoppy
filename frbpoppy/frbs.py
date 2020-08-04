"""Class to hold FRB source properties."""
import numpy as np
import pandas as pd


class FRBs:
    """Class containing FRB properties."""

    def __init__(self):
        """Initializing."""
        # Location properties
        self.ra = None  # Right ascension [frac degree]
        self.dec = None  # Declination [frac degree]
        self.dist_co = None  # Comoving distance [Gpc]
        self.gb = None  # Galactic latitude [frac degree]
        self.gl = None  # Galactic longitude [frac degree]
        self.gx = None  # Galactic X coordinate [Gpc]
        self.gy = None  # Galactic Y coordinate [Gpc]
        self.gz = None  # Galactic Z coordinate [Gpc]
        self.z = 0  # Redshift

        # Dispersion measure properties
        self.dm_host = 0  # DM host galaxy [pc/cm^3]
        self.dm_igm = 0  # DM intergalactic medium [pc/cm^3]
        self.dm_mw = 0  # DM Milky Way [pc/cm^3]
        self.dm = None  # Total DM [pc/cm^3]

        # Intrinsic properties
        self.lum_bol = None  # Isotropic equivalent bolometric lum. [ergs/s]
        self.si = None  # Spectral index
        self.w_arr = None  # Pulse width at Earth [ms]
        self.w_int = None  # Intrinsic pulse width [ms]

        # Repeat properties
        self.time = None  # Burst time stamps [days]

        # Detection properties
        self.fluence = None  # Fluence [Jy ms]
        self.offset = None  # Offset from beam centre [frac deg]
        self.s_peak = None  # Peak flux density [Jy]
        self.snr = None  # Signal to Noise ratio
        self.t_dm = 0  # Dispersion meausre smearing timescale [ms]
        self.t_scat = 0  # Scattering timescale [ms]
        self.T_sky = 0  # Sky temperature [K]
        self.T_sys = 0  # Total system temperature [K]
        self.w_eff = None  # Effective pulse width [ms]

        # Software properties
        self.index = None  # Index to keep track of FRBs

    def __str__(self):
        """Define how to print an FRB object to a console."""
        s = 'FRBs properties:'
        for key, value in vars(self).items():
            if isinstance(value, np.ndarray):
                value = f'{len(value)} elements - {value[:2]} etc.'
            s = '\n\t'.join([s, f"{key}: {value}"])

        return s

    def apply(self, mask):
        """Apply a Numpy array to all parameters.

        Args:
            mask (array): Masking array to apply to all frb parameters.

        """
        for attr in self.__dict__.keys():
            parm = getattr(self, attr)
            if type(parm) is np.ndarray:
                # 1D mask on 1D or 2D array
                if mask.ndim == 1:
                    setattr(self, attr, parm[mask])
                else:
                    # 2D mask on 1D array
                    if parm.ndim == 1:
                        setattr(self, attr, parm[mask.any(axis=1)])
                    # 2D mask on 2D array
                    else:
                        # Apply mask in standard way
                        parm_nan = np.where(mask, parm, np.nan)
                        # Keep rows with any True values
                        row_mask = np.any(mask, axis=1)
                        parm_nan = parm_nan[row_mask, :]
                        # Keep columns with any True values
                        col_mask = np.any(mask[row_mask], axis=0)
                        parm_nan = parm_nan[:, col_mask]
                        # Set attribute
                        setattr(self, attr, parm_nan)

    def clean_up(self):
        """Clean up 2D parameter arrays by left justifying them."""
        # First apply a time mask everywhere
        if type(self.time) is np.ndarray:
            time_mask = ~np.isnan(self.time)
            self.apply(time_mask)

        # Then left justify all 2D parameters
        for attr in self.__dict__.keys():
            parm = getattr(self, attr)
            if (type(parm) is np.ndarray) and (parm.ndim == 2):
                mask = ~np.isnan(parm)
                justified_mask = np.sort(mask, 1)[:, ::-1]
                out = np.full(parm.shape, np.nan)
                out[justified_mask] = parm[mask]
                out = out[:, ~np.all(np.isnan(out), axis=0)]
                setattr(self, attr, out)

    def to_df(self):
        """Convert properties to a Pandas DataFrame."""
        # Find all source properties
        df = pd.DataFrame()
        for attr in self.__dict__.keys():
            parm = getattr(self, attr)
            if type(parm) is not np.ndarray:
                continue

            # 2D arrays to 1D
            if parm.ndim > 1:
                if parm.size > 0:
                    df[attr] = parm[~np.isnan(self.time)]
            # 1D array to match length 2D->1D arrays
            elif type(self.time) is np.ndarray:
                concate = np.array([parm, ]*self.time.shape[1]).transpose()
                if concate.size > 0:
                    df[attr] = concate[~np.isnan(self.time)]
            else:
                df[attr] = parm

        return df

    def to_csv(self, path):
        """Export properties to a csv file."""
        self.to_df().to_csv(path)

    def calc_logn_logs(self, parameter='fluence', min_p=None, max_p=None):
        """Rough method to retrive the slope of a LogNlogS distribution."""
        p = getattr(self, parameter)

        # For repeaters with 2D arrays of parameters
        if p.ndim > 1:
            p = p.flatten()

        # Add value limits if wished
        if min_p is None:
            f_0 = np.min(p)
        else:
            f_0 = min_p
            p = p[p >= min_p]

        if max_p is not None:
            p = p[p <= max_p]

        n = p.size
        alpha = -1/((1/n)*sum([np.log(f/f_0) for f in p]))
        alpha *= (n-1)/n  # Removing bias in alpha
        alpha_err = n*alpha/((n-1)*(n-2)**0.5)
        norm = n / (f_0**alpha)  # Normalisation at lowest value

        return alpha, alpha_err, norm


if __name__ == '__main__':
    frbs = FRBs()
    import IPython; IPython.embed()
