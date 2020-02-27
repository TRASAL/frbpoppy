"""Spectral index distributions."""
import numpy as np
import frbpoppy.gen_dists as gd


def constant(value=1e40, shape=1):
    """Good for standard candles."""
    return np.full(shape, value)


def powerlaw(low=1e40, high=1e45, power=0, shape=1):
    """Draw luminosities from powerlaw distribution."""
    return gd.powerlaw(low=low, high=high, power=power, shape=shape)


def gauss(mean=1e35, std=1e2, shape=1):
    """Generate luminosties from a Gaussian/Normal distribution.

    Args:
        mean (float): Mean luminosity [ergs/s]
        std (float): Standard deviation [ergs/s]
        shape (tuple): Required array shape

    Returns:
        array: Luminosties

    """
    return gd.trunc_norm(mean, std, shape)


def log10normal(mean=1e40, std=1e2, shape=1):
    """Draw luminosity from a log10normal distribution."""
    return gd.log10normal(np.log10(mean), np.log10(std), shape)
