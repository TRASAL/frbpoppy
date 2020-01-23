"""Spectral index distributions."""
import numpy as np


def gauss(mu=-1.4, sigma=1, shape=1):
    """Generate spectral indices from a Gaussian distribution.

    Args:
        mu (float): Mean spectral index
        sigma (float): Spread of spectral index
        shape (tuple): Required array shape

    Returns:
        array: spectral indices

    """
    return np.random.normal(mu, sigma, shape).astype(np.float32)


def gauss_per_source(b_sigma=0.05, dist=gauss, shape=(1, 1), z=0, **kwargs):
    """Distribute spectral indices per source using a Gaussian function.

    Generate bursts per source using a given distribution, then use those as
    the mean for the distribution for more bursts from that source
    """
    mu = dist(shape=shape[0], **kwargs)

    # Check whether multiple bursts per source necessary
    if len(shape) < 2:
        return mu

    shape = (shape[1], shape[0])
    return np.random.normal(mu, b_sigma, shape).astype(np.float32).T
