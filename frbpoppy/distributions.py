"""Define distributions from which to get random numbers."""
import numpy as np
import math
from scipy.stats import truncnorm
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
        return 10**np.random.uniform(math.log10(low), math.log10(high), n_gen)

    def sample(n_gen):
        pl = np.random.uniform(0, 1, n_gen)**(1/power)
        if power > 0:
            addition = math.log10(high)
        else:
            addition = math.log10(low)
        log_pl = np.log10(pl) + addition
        pl = 10**log_pl
        return pl

    def accept(pl):
        if power > 0:
            return pl >= low
        else:
            return pl <= high

    pl = sample(n_gen)
    mask = accept(pl)
    reject, = np.where(~mask)
    while reject.size > 0:
        fill = sample(reject.size)
        mask = accept(fill)
        pl[reject[mask]] = fill[mask]
        reject = reject[~mask]

    return pl


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


def trunc_norm(mu, sigma, n_gen=1, lower=0, upper=np.inf):
    """Draw from a truncated normal distribution.

    Args:
        mu (number): Mu
        sigma (number): Sigma
        n_gen (number): Number to generate
        lower (number): Lower limit
        upper (number): Upper limit

    Returns:
        array: Numpy of required length

    """
    if sigma == 0:
        return np.full(n_gen, mu)
    left = (lower-mu)/sigma
    right = (upper-mu)/sigma
    d = truncnorm.rvs(left, right, loc=mu, scale=sigma, size=n_gen)
    return d
