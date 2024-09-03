"""Functions to generate dispersion measure distributions."""
import numpy as np
#import numexpr as ne
import frbpoppy.gen_dists as gd

#global rng
#rng = np.random.default_rng()

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
            ('normal', 'lognormal', 'log10normal')

    Returns:
        dm_igm (array): Dispersion measure of intergalactic medium [pc/cm^3]

    """
    if std is None:
        std = 0.2*slope*z
        #std = ne.evaluate("0.2*slope*z")

    # Set up spread distribution
    mean = slope*z
    #mean = ne.evaluate("slope*z")
    if spread_dist == 'normal':
        f = np.random.normal
        #f = rng.normal
    elif spread_dist == 'lognormal':
        def f(mean, std):
            return gd.lognormal(mean, std)
    elif spread_dist == 'log10normal':
        def f(mean, std):
            return gd.log10normal(mean, std)
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
    #return ne.evaluate("dm_host / (1 + z)")


def log10normal(mean=100, std=200, n_srcs=1, z=0):
    """Generate a log10 normal dm host distribution.

    Args:
        mean (float): Mean DM [pc/cm^3].
        std (float): Standard deviation DM [pc/cm^3].
        n_srcs (int): Number of sources for which to generate values.
        z (int): Redshift of sources.

    Returns:
        array: DM host [pc/cm^3]

    """
    dm_host = gd.log10normal(mean, std, n_srcs).astype(np.float32)
    return dm_host / (1 + z)
    #return ne.evaluate("dm_host / (1 + z)")


def lognormal(mean=100, std=200, n_srcs=1, z=0):
    """Generate a lognormal dm host distribution.

    Args:
        mean (float): Mean DM [pc/cm^3].
        std (float): Standard deviation DM [pc/cm^3].
        n_srcs (int): Number of sources for which to generate values.
        z (int): Redshift of sources.

    Returns:
        array: DM host [pc/cm^3]

    """
    dm_host = gd.lognormal(mean, std, n_srcs).astype(np.float32)
    return dm_host / (1 + z)
    #return ne.evaluate("dm_host / (1 + z)")
