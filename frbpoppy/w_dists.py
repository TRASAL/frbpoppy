"""Pulse width distributions."""
import numpy as np

import frbpoppy.gen_dists as gd


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
    # These are methods to deal with numpy dimensionality aspects
    elif np.argmax(w_int.shape) != np.argmax(z.shape):
        return w_int*(1+z[:, None]).T
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


def gauss(mean=1, std=2, shape=1, z=0):
    """Generate pulse widths from a Gaussian/Normal distribution.

    Args:
        mean (float): Mean pulse width [ms]
        std (float): Standard deviation of the pulse widths [ms]
        shape (tuple): Required array shape
        z (array): Redshift of pulses

    Returns:
        tuple: intrinsic pulse widths, pulse widths at Earth

    """
    w_int = gd.trunc_norm(mean, std, shape).astype(np.float32)
    w_arr = calc_w_arr(w_int, z=z)
    return w_int, w_arr


def log10normal(mean=0.1, std=0.5, shape=1, z=0):
    """Draw burst from log10normal distribution.

    Args:
        shape (tuple): Required array shape
        z (array): Redshift of pulses

    Returns:
        type: Description of returned object.

    """
    w_int = gd.log10normal(mean, std, shape).astype(np.float32)
    w_arr = calc_w_arr(w_int, z=z)
    return w_int, w_arr
