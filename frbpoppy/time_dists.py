"""Distributions with which burst times can be generated."""
import numpy as np
from scipy.special import gamma


def single(n_srcs=1, n_days=1, z=0):
    """Generate a series of one-off burst times.

    Args:
        n_srcs (int): Number of sources
        n_days (int): Number of days
        z (array): Redshift of sources
    """
    time = np.random.uniform(0, n_days, n_srcs).astype(np.float32)
    time *= (1+z)
    return time[:, np.newaxis]


def regular(rate=2, n_srcs=1, n_days=1, z=0):
    """Generate a series of regular spaced burst times.

    Args:
        rate (float): Number of events per day
        n_srcs (int): Number of sources
        n_days (int): Number of days
        z (array): Redshift of sources
    """
    time_range = np.arange(0, n_days, step=1/rate, dtype=np.float32)

    # Copy to multiple sources
    time = np.tile(time_range, (n_srcs, 1))

    # Add redshift
    time *= np.atleast_1d(1+z)[:, np.newaxis]
    time[(time > n_days)] = np.nan

    return time


def cyclic(rate=2, n_days=1, n_srcs=1, period=1, frac=.1, z=0):
    """Generate a series of uniform burst times within an activity cycle."

    Args:
        rate (float/array): Number of events per day
        n_days (int): Number of days
        n_srcs (int): Number of sources
        period (float/array): Period of activity cycle (days)
        frac (float/array): Fraction of activity cycle a source is active
        z (float/array): Redshift of sources
    """
    # ensure arguments that may be floats or arrays are arrays (length is number of sources)
    rate = np.atleast_1d(rate)
    period = np.atleast_1d(period)
    frac = np.atleast_1d(frac)
    z = np.atleast_1d(z)

    # get number of bursts to generate for each source
    nburst_per_cycle = (rate * frac * period).astype(int)
    nburst_per_source = ((n_days / period) * nburst_per_cycle).astype(int)

    # generate burst arrival times within each active period
    times = np.random.uniform(0, frac[:, np.newaxis] * period[:, np.newaxis],
                              (n_srcs, nburst_per_source.max())).astype(np.float32)

    # need to add jump when going to next cycle
    for src_ind in range(n_srcs):

        times += (np.arange(nburst_per_source.max()) // nburst_per_cycle[src_ind]) * period[src_ind]

    # Add redshift
    t0 = times[:, 0][:, np.newaxis]
    times -= t0
    times *= (1+z)[:, np.newaxis]
    times += t0

    # Remove bursts past n_days
    times[(times > n_days)] = np.nan

    return times


def poisson(**kwargs):
    """Generate a series of poisson times.

    Wrapper to iteratively generate from _poisson_dist

    Args:
        Kwargs from iteratively_gen_times
    """
    # Set distribution from which to draw
    return iteratively_gen_times(_poisson_dist, **kwargs)


def _poisson_dist(dims, rate=0.1):
    """Draw values from a poissonian distribution.

    Args:
        rate (float): Expected number of events per day.
    """
    if not isinstance(rate, np.ndarray):
        return np.random.exponential(1/rate, dims).astype(np.float32)
    else:  # Allow for an array of lambdas
        dims = dims[::-1]
        return np.random.exponential(1/rate, dims).astype(np.float32).T


def clustered(**kwargs):
    """Generate burst times following a Weibull distribution.

    Wrapper to iteratively generate from _weibull_dist

    Args:
        Kwargs from iteratively_gen_times
    """
    return iteratively_gen_times(_weibull_dist, **kwargs)


def _weibull_dist(dims, r=5.7, k=0.34):
    """Generate burst times following a Weibull distribution.

    Args:
        r (float): Rate parameter
        k (float): Shape parameter
    """
    lam = 1/(r*gamma(1 + 1/k))
    if not any([isinstance(p, np.ndarray) for p in (r, k)]):
        return lam*np.random.weibull(k, dims).astype(np.float32)
    else:  # Allow for an array for r's
        dims = dims[::-1]
        return (lam*np.random.weibull(k, dims).astype(np.float32)).T


def iteratively_gen_times(dist, n_srcs=1, n_days=1, z=0, **kwargs):
    """Generate burst times in an iterative manner using a distribution.

    Args:
        dist (func): Distribution from which to draw burst times
        n_srcs (int): Number of sources
        n_days (int): Number of days
        z (array): Redshift of sources
    """
    # Determine the maximum possible number of bursts per source to include
    log_size = 1
    m = int(10**log_size)
    dims = (n_srcs, m)
    time = dist(dims, **kwargs)

    # Converting intervals to time stamps
    time = np.cumsum(time, axis=1)  # This is in fraction of days

    if isinstance(z, np.ndarray):
        time = time*(1+z)[:, np.newaxis]
    else:
        time = time*(1+z)

    # Mask any frbs over the maximum time (Earth perspective)
    time[(time > n_days)] = np.nan

    # Iteratively add extra bursts until over limit
    mask = ~np.isnan(time[:, -1])  # Where more bursts are needed
    sum_mask = np.count_nonzero(mask)

    while sum_mask != 0:  # Add additional columns
        m = int(10**log_size)
        ms = np.full((n_srcs, m), np.nan)

        # Ensure the size of the keyword values for dist are correct
        new_kwargs = kwargs.copy()
        for kw, value in kwargs.items():
            if isinstance(value, np.ndarray) and sum_mask != dims[0]:
                new_kwargs[kw] = kwargs[kw][mask]

        new = dist((sum_mask, m), **new_kwargs)
        new = np.cumsum(new, axis=1)

        # Add redshift correction
        masked_z = z
        if isinstance(z, np.ndarray):
            masked_z = z[mask][:, np.newaxis]
        new *= (1+masked_z)

        new += time[:, -1][mask][:, np.newaxis]  # Ensure cumulative
        new[(new > n_days)] = np.nan  # Apply filter
        ms[mask] = new  # Set up additional columns
        time = np.hstack((time, ms))  # Add to original array

        # Set up for next loop
        if sum_mask == n_srcs:
            log_size += 0.5
        mask = ~np.isnan(time[:, -1])
        sum_mask = np.count_nonzero(mask)

    # Remove any columns filled with NaNs
    time = time[:, ~np.all(np.isnan(time), axis=0)]

    return time
