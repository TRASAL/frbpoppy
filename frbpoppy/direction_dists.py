"""Spatial direction distributions of frb sources."""
import numpy as np


def uniform(min_ra=0, max_ra=360, min_dec=-90, max_dec=90, n_srcs=1):
    # Add random directional coordinates
    u = np.random.uniform
    ra = u(min_ra, max_ra, n_srcs)
    dec = np.rad2deg(np.arccos(u(min_dec/90, max_dec/90, n_srcs))) - 90
    return ra, dec
