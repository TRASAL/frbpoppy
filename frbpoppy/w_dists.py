"""Pulse width distributions."""
import numpy as np


def calc_w_arr(w_int, z=0):
    """Calculate the pulse width upon arrival to Earth.

    Args:
        w_int (array): Intrinsic pulse duration
        z (array): Redshift of pulses

    Returns:
        array: Pulse width at arrival to Earth

    """
    if w_int.ndim == 1 or isinstance(z, int) or isinstance(z, float):
        return w_int*(1+z)
    else:
        return w_int*(1+z[:, None])


def constant(value=1., shape=1, z=0):
    """Generate pulse widths at a constant value."""
    w_int = np.full(shape, value).astype(np.float32)
    return w_int, calc_w_arr(w_int, z=z)


def uniform(low=0, high=10, shape=1, z=0):
    """Generate pulse widths from a uniform distribution.

    Args:
        low (float): Minimum pulse width [ms]
        high (float): Maximum pulse width [ms]
        shape (tuple): Required array shape
        z (array): Redshift of pulses

    Returns:
        type: Description of returned object.

    """
    w_int = np.random.uniform(low, high, shape).astype(np.float32)
    w_arr = calc_w_arr(w_int, z=z)
    return w_int, w_arr


def lognormal(mu=0.1, sigma=0.5, shape=1, z=0):
    """Short summary.

    Args:
        shape (tuple): Required array shape
        z (array): Redshift of pulses

    Returns:
        type: Description of returned object.

    """
    w_int = np.random.lognormal(mu, sigma, shape).astype(np.float32)
    w_arr = calc_w_arr(w_int, z=z)
    return w_int, w_arr


def gauss_per_source(src_sigma=0.05, dist=uniform, shape=(1, 1), z=0,
                     **kwargs):
    """Distribute pulse widths per source using a Gaussian function.

    Generate bursts per source using a given distribution, then use those as
    the mean for the distribution for more bursts from that source
    """
    mu, s = dist(shape=shape[0], z=z, **kwargs)

    # Check whether multiple bursts per source needed
    if len(shape) < 2:
        return mu, s

    shape = (shape[1], shape[0])
    w_int = np.random.normal(mu, src_sigma, shape).astype(np.float32).T
    w_arr = calc_w_arr(w_int, z=z)
    return w_int, w_arr
