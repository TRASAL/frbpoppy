"""Spectral index distributions."""
import numpy as np
import frbpoppy.gen_dists as gd


def powerlaw(low=1e40, high=1e45, power=0, shape=1):
    """Draw luminosities from powerlaw distribution."""
    return gd.powerlaw(low=low, high=high, power=power, shape=shape)


def gauss_per_source(src_sigma=0.05, dist=powerlaw, shape=(1, 1), z=0,
                     **kwargs):
    """Distribute spectral indices per source using a Gaussian function.

    Generate bursts per source using a given distribution, then use those as
    the mean for the distribution for more bursts from that source

    Args:
        src_sigma: Standard deviation over bursts.
        dist: Distribution from which to draw source properties.
        shape: Wanted shape of array of spectral indices.
        z: Redshift.
    """
    mu = dist(shape=shape[0], **kwargs)

    # Check whether multiple bursts per source necessary
    if len(shape) < 2:
        return mu

    shape = (shape[1], shape[0])

    return np.random.normal(mu, src_sigma, shape).astype(np.float32).T
