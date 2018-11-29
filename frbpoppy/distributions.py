"""Define distributions from which to get random numbers."""
import numpy as np
import random
import frbpoppy.precalc as pc


def powerlaw(low, high, power, n_gen=1):
    """
    Return random variables distributed according to power law.

    The power law distribution power is simply the power, not including a minus
    sign (P scales with x^n with n the power). A flat powerlaw can therefore
    be created by taking setting power to zero.

    Args:
        low (float): Lower limit of distribution
        high (float): Higher limit of distribution
        power (float): Power of power law distribution
        n_gen (int): Number of values to be generated

    Returns:
        array: Random variable picked from power law distribution

    """
    if low > high:
        low, high = high, low

    if power == 0:
        return 10**np.random.uniform(np.log10(low), np.log10(high), n_gen)

    return np.random.uniform(low, high, n_gen)**(1/power)


def z_from_sfr(z_max=2.5, n_gen=1):
    """
    Return a random redshift for sources following the Star Formation Rate.

    Follows Madau & Dickinson (2014), eq. 15. For more info see
    https://arxiv.org/pdf/1403.0007.pdf

    """
    def sfr(z):
        return (1+z)**2.7/(1+((1+z)/2.9)**5.6)

    def sample(n_gen):
        return np.random.uniform(0, z_max, (n_gen,))

    def accept(x):
        return np.random.rand(x.size)*9.0 <= sfr(x)

    z = sample(n_gen)
    mask = accept(z)
    reject, = np.where(~mask)
    while reject.size > 0:
        fill = sample(reject.size)
        mask = accept(fill)
        z[reject[mask]] = fill[mask]
        reject = reject[~mask]

    return z


def z_from_csmd(z_max=6.0, n_gen=1):
    """
    Return a random redshift for sources following the Stellar Mass Density.

    Follows Madau & Dickinson (2014), eq. 2 & 15. For more info see
    https://arxiv.org/pdf/1403.0007.pdf

    """
    # z = None
    #
    # while not z:
    #     x = random.random()*z_max
    #     y = random.random()*0.00065
    #     if y <= pc.csmd_table(x):
    #         z = x
    csmd = pc.CSMDTable().lookup

    def sample(n_gen):
        return np.random.uniform(0, z_max, (n_gen,))

    def accept(x):
        return np.random.rand(x.size)*0.00065 <= csmd(x)

    # Accept - rejct algorithm
    z = sample(n_gen)
    mask = accept(z)
    reject, = np.where(~mask)
    while reject.size > 0:
        fill = sample(reject.size)
        mask = accept(fill)
        z[reject[mask]] = fill[mask]
        reject = reject[~mask]

    return z
