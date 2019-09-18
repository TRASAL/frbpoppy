# -*- coding: future_fstrings -*-
"""Define distributions from which to get random numbers."""
import numpy as np
import math
from scipy.stats import truncnorm
import frbpoppy.precalc as pc


def powerlaw(low, high, power, shape=1):
    """
    Return random variables distributed according to power law.

    The power law distribution power is simply the power, not including a minus
    sign (P scales with x^n with n the power). A flat powerlaw can therefore
    be created by taking setting power to zero.

    Args:
        low (float): Lower limit of distribution
        high (float): Higher limit of distribution
        power (float): Power of power law distribution
        shape (int/tuple): Shape of array to be generated. Can also be a int.

    Returns:
        array: Random variable picked from power law distribution

    """
    if low > high:
        low, high = high, low

    if power == 0 or low == high:
        return 10**np.random.uniform(np.log10(low), np.log10(high), shape)

    def sample(n_gen):
        pl = np.random.uniform(0, 1, n_gen)**(1/power)
        if power > 0:
            addition = np.log10(high)
        else:
            addition = np.log10(low)
        log_pl = np.log10(pl) + addition
        pl = 10**log_pl
        return pl

    def accept(pl):
        if power > 0:
            return pl >= low
        else:
            return pl <= high

    if isinstance(shape, tuple):
        n_gen = shape[0] * shape[1]
    else:
        n_gen = shape

    pl = sample(n_gen)
    mask = accept(pl)
    reject, = np.where(~mask)
    while reject.size > 0:
        fill = sample(reject.size)
        mask = accept(fill)
        pl[reject[mask]] = fill[mask]
        reject = reject[~mask]

    if isinstance(shape, tuple):
        pl = pl.reshape(shape)
        
    return pl


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


def oppermann_pen():
    """Following Oppermann & Pen (2017), simulate repeat times.
    Returns:
        ts (list): List of burst times
    """
    # TODO Numpy-ify
    r = 5.7
    k = 0.34

    ts = []
    t_tot = 0.5  # Assuming a maximum of 12 hours on one spot
    t_sum = 0.0

    # Get time of bursts
    while t_sum < t_tot:
        t = r*np.random.weibull(k)
        t_sum += t
        ts.append(t_sum)

    # Convert to seconds
    ts = [t*86400 for t in ts[:-1]]

    return ts
