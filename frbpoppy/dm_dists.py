"""Functions to generate Dispersion measure distributions."""
import numpy as np
import frbpoppy.gen_dists as gd


def constant(value=100, n_srcs=1):
    """Adopt a constant DM value similar to Thorton."""
    return np.full(n_srcs, value).astype(np.float32)


def ioka(z=0, slope=950, sigma=None, spread_dist='normal'):
    """Calculate the contribution of the igm to the dispersion measure.

    Follows Ioka (2003) and Inoue (2004), with default slope value falling
    in between the Cordes and Petroff reviews.

    Args:
        z (array): Redshifts.
        slope (float): Slope of the DM-z relationship.
        sigma (float): Spread around the DM-z relationship.
        spread_dist (str): Spread function option. Choice from
            ('normal', 'lognormal')

    Returns:
        dm_igm (array): Dispersion measure of intergalactic medium [pc/cm^3]

    """
    if sigma is None:
        sigma = 0.2*slope*z

    # Set up spread distribution
    if spread_dist == 'normal':
        f = np.random.normal
    elif spread_dist == 'lognormal':
        f = np.random.lognormal
    else:
        raise ValueError('spread_dist input not recognised')
    return f(slope*z, sigma).astype(np.float32)


def gauss(mu=100, sigma=200, n_srcs=1, z=0):
    """Generate dm host contributions similar to Tendulkar.

    Args:
        mu (float): Mean DM [pc/cm^3].
        sigma (float): Standard deviation DM [pc/cm^3].
        n_srcs (int): Number of sources for which to generate values.
        z (int): Redshift of sources.

    Returns:
        array: DM host [pc/cm^3]

    """
    dm_host = gd.trunc_norm(mu, sigma, n_srcs).astype(np.float32)
    return dm_host / (1 + z)


def lognormal(mu=100, sigma=200, n_srcs=1, z=0):
    """Generate a log normal dm host distribution.

    Args:
        mu (float): Mean DM [pc/cm^3].
        sigma (float): Standard deviation DM [pc/cm^3].
        n_srcs (int): Number of sources for which to generate values.
        z (int): Redshift of sources.

    Returns:
        array: DM host [pc/cm^3]

    """
    dm_host = np.random.lognormal(mu, sigma, n_srcs).astype(np.float32)
    return dm_host / (1 + z)
