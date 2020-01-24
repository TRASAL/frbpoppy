"""Spatial direction distributions for FRB sources."""
import numpy as np


def uniform(min_ra=0, max_ra=360, min_dec=-90, max_dec=90, n_srcs=1):
    """Generate a uniform distribution of pointings.

    Args:
        min_ra (float): Minimum right ascenion [frac deg].
        max_ra (float): Maximum right ascenion [frac deg].
        min_dec (float): Minimum declination [frac deg].
        max_dec (float): Maximum declination [frac deg].
        n_srcs (int): Number of sources for which to generate.

    Returns:
        tuple: RA, Dec arrays [frac deg]

    """
    u = np.random.uniform
    ra = u(min_ra, max_ra, n_srcs)
    dec = np.rad2deg(np.arccos(u(min_dec/90, max_dec/90, n_srcs))) - 90
    return ra, dec
