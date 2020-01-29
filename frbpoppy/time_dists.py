"""Distributions with which burst times can be generated."""
import numpy as np
from scipy.special import gamma


def regular(lam=2, n_srcs=1, n_days=1, z=0):
    """Generate a series of regular spaced burst times.

    Args:
        lam (float): Number of events per day
        n_srcs (int): Number of sources
        n_days (int): Number of days
        z (array): Redshift of sources
    """
    time_range = np.linspace(0, n_days, n_days*lam, endpoint=False,
                             dtype=np.float32)

    # Copy to multiple sources
    time = np.tile(time_range, (n_srcs, 1))

    # Add random offsets
    time_offset = np.random.uniform(0, time[:, 1])
    time += time_offset[:, np.newaxis]

    # Add redshift
    time = time*(1+z)[:, np.newaxis]
    time[(time > n_days)] = np.nan

    return time


def poisson(lam=0.1, **kwargs):
    """Generate a series of poisson times.

    Args:
        lam (float): Expected number of events per day.
        Kwargs from iteratively_gen_times
    """
    # Set distribution from which to draw
    def poisson(dims, lam=lam):
        return np.random.exponential(1/lam, dims).astype(np.float32)

    return iteratively_gen_times(poisson, **kwargs)


def clustered(r=5.7, k=0.34, **kwargs):
    """Generate burst times following a Weibull distribution.

    Args:
        r (float): Rate parameter
        k (float): Shape parameter
        Kwargs from iteratively_gen_times
    """
    # Set distribution from which to draw
    def weibull(dims, r=r, k=k):
        lam = 1/(r*gamma(1 + 1/k))
        return lam*np.random.weibull(k, dims).astype(np.float32)

    return iteratively_gen_times(weibull, **kwargs)


def iteratively_gen_times(dist, n_srcs=1, n_days=1, z=0):
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
    time = dist(dims)
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
        new = dist((sum_mask, m))
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
