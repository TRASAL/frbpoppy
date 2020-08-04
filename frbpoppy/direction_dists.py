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
    min_dec_x = np.cos(np.deg2rad(min_dec + 90))
    max_dec_x = np.cos(np.deg2rad(max_dec + 90))
    dec = np.rad2deg(np.arccos(u(max_dec_x, min_dec_x, n_srcs))) - 90
    return ra, dec


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    ra, dec = uniform(min_ra=0, max_ra=360, min_dec=-45, max_dec=90,
                      n_srcs=int(1e4))
    ra[ra > 180] = ra[ra > 180] - 360  # Wrap for plotting
    plt.subplot(111, projection="aitoff")
    plt.scatter(np.radians(ra), np.radians(dec))
    plt.grid(True)
    plt.show()
