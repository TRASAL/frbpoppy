"""Spectral index distributions."""
import numpy as np


def constant(value=-1.4, shape=1):
    """Good for adopting a single value."""
    return np.full(shape, value).astype(np.float32)


def gauss(mean=-1.4, std=1, shape=1):
    """Generate spectral indices from a Gaussian distribution.

    Args:
        mean (float): Mean spectral index
        std (float): Spread of spectral index
        shape (tuple): Required array shape

    Returns:
        array: spectral indices

    """
    return np.random.normal(mean, std, shape).astype(np.float32)


def gauss_per_source(b_std=0.05, dist=gauss, shape=(1, 1), z=0, **kwargs):
    """Distribute spectral indices per source using a Gaussian function.

    Generate bursts per source using a given distribution, then use those as
    the mean for the distribution for more bursts from that source
    """
    mean = dist(shape=shape[0], **kwargs)

    # Check whether multiple bursts per source necessary
    if len(shape) < 2:
        return mean

    shape = (shape[1], shape[0])
    return np.random.normal(mean, b_std, shape).astype(np.float32).T
