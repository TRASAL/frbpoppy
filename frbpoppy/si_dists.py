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
