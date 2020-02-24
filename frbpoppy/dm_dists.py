"""Functions to generate Dispersion measure distributions."""
import numpy as np
import frbpoppy.gen_dists as gd
from frbpoppy.misc import lognormal_to_normal


def constant(value=100, n_srcs=1):
    """Adopt a constant DM value similar to Thorton."""
    return np.full(n_srcs, value).astype(np.float32)


def ioka(z=0, slope=950, std=None, spread_dist='normal'):
    """Calculate the contribution of the igm to the dispersion measure.

    Follows Ioka (2003) and Inoue (2004), with default slope value falling
    in between the Cordes and Petroff reviews.

    Args:
        z (array): Redshifts.
        slope (float): Slope of the DM-z relationship.
        std (float): Spread around the DM-z relationship.
        spread_dist (str): Spread function option. Choice from
            ('normal', 'lognormal')

    Returns:
        dm_igm (array): Dispersion measure of intergalactic medium [pc/cm^3]

    """
    if std is None:
        std = 0.2*slope*z

    # Set up spread distribution
    mean = slope*z
    if spread_dist == 'normal':
        f = np.random.normal
    elif spread_dist == 'lognormal':
        f = np.random.lognormal
        mean, std = lognormal_to_normal(mean, std)
    else:
        raise ValueError('spread_dist input not recognised')
    return f(mean, std).astype(np.float32)


def gauss(mean=100, std=200, n_srcs=1, z=0):
    """Generate dm host contributions similar to Tendulkar.

    Args:
        mean (float): Mean DM [pc/cm^3].
        std (float): Standard deviation DM [pc/cm^3].
        n_srcs (int): Number of sources for which to generate values.
        z (int): Redshift of sources.

    Returns:
        array: DM host [pc/cm^3]

    """
    dm_host = gd.trunc_norm(mean, std, n_srcs).astype(np.float32)
    return dm_host / (1 + z)


def lognormal(mean=100, std=200, n_srcs=1, z=0):
    """Generate a log normal dm host distribution.

    Args:
        mean (float): Mean DM [pc/cm^3].
        std (float): Standard deviation DM [pc/cm^3].
        n_srcs (int): Number of sources for which to generate values.
        z (int): Redshift of sources.

    Returns:
        array: DM host [pc/cm^3]

    """
    mean, std = lognormal_to_normal(mean, std)
    dm_host = np.random.lognormal(mean, std, n_srcs).astype(np.float32)
    return dm_host / (1 + z)
