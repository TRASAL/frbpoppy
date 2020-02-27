"""Distributions with which burst times can be generated."""
import numpy as np
from tqdm import tqdm
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

    # Add random offsets
    if time.shape[1] > 1:
        time_offset = np.random.uniform(0, time[:, 1])
    else:
        time_offset = np.random.uniform(0, time[:, 0])
    time += time_offset[:, np.newaxis]

    # Add redshift
    time = time*(1+z)[:, np.newaxis]
    time[(time > n_days)] = np.nan

    return time


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

    time = np.cumsum(time, axis=1)  # This is in fraction of days
    # The first burst time is actually the time since the previous one
    # You want to be at in a random time in between those
    time_offset = np.random.uniform(0, time[:, 0])
    time -= time_offset[:, np.newaxis]

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
