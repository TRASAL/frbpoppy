"""Define general distributions from which to get random numbers."""
import numpy as np
from scipy.stats import truncnorm


def powerlaw(low, high, power, shape=1):
    """
    Return random variables distributed according to power law.

    The power law distribution power is simply the power, not including a minus
    sign (N(x) scales with x^p with p the power). A flat powerlaw can therefore
    be created by taking setting power to zero.

    Args:
        low (float): low limit of distribution
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


def trunc_norm(mean, std, shape=1, low=0, high=np.inf):
    """Draw from a truncated normal distribution.

    Args:
        mean (number): Mean of the normal distribution.
        std (number): Standard deviation of the distribution.
        shape (number): Number to generate.
        low (number): Lower limit.
        high (number): Higher limit.

    Returns:
        array: Numpy of required length

    """
    if not isinstance(std, np.ndarray) and std == 0:
        return np.full(shape, mean)
    left = (low-mean)/std
    right = (high-mean)/std
    return truncnorm.rvs(left, right, loc=mean, scale=std, size=shape)


def log10normal(mean, std, shape):
    """Random values from a normal distribution in the log space.

    I could never quite figure out what to expect of a lognormal distribution,
    so I implemented a log10normal distribution which looks like a normal
    distribution when binned in the log10 space.
    """
    mean, std = np.log10(mean), np.log10(std)
    return 10**np.random.normal(mean, std, shape)


def calc_lognormal_input(mean_x, std_x):
    """Calculate the mean and std of a lognormal distribution.

    See
    https://en.wikipedia.org/wiki/Log-normal_distribution
    """
    normal_std = np.sqrt(np.log(1 + (std_x**2/mean_x)**2))
    normal_mean = np.log(mean_x**2 / np.sqrt(mean_x**2 + std_x**2))
    return normal_mean, normal_std


def lognormal(mean, std, shape):
    """Calculate the mean and std from the underlying distribution.

    See
    https://en.wikipedia.org/wiki/Log-normal_distribution
    """
    mean, std = calc_lognormal_input(mean, std)
    return np.random.lognormal(mean, std, shape)
