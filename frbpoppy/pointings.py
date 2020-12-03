"""Functions to generate pointings."""
import numpy as np
from datetime import datetime, timedelta

import frbpoppy.galacticops as go


def transit(n_gen, lat=None, lon=None, t_obs=None):
    """Generate RA and Dec pointing coordinates for a transit telescope.

    Args:
        n_gen (int): Number of pointings wanted.
        lat (float): Latitude of telescope.
        lon (float): Longitude of telescope (minus for west).
        t_obs (type): Time for each pointing.

    Returns:
        ra, dec: Numpy arrays with coordinates of pointings

    """
    # Pointings are randomly placed between the year 2000 and 2100
    date_min = go.random_date(datetime(2000, 1, 1), datetime(2100, 1, 1))
    date_max = date_min + timedelta(seconds=int(t_obs*n_gen))
    time_delta = np.timedelta64(int(t_obs), 's')
    times = np.arange(date_min, date_max, time_delta, dtype='datetime64')
    ra = (go.datetime_to_gmst(times) + lon) % 360
    dec = np.ones(n_gen)*lat

    return ra, dec


def tracking(n_gen, **kwargs):
    """Generate RA, Dec pointing coordinates for a tracking telescope.

    Uses a try-accept algorithm to generate pointings with a survey. For
    details on the sunflower algoritm used to distribute pointings see
    See https://stackoverflow.com/a/44164075/11922471. Pointings are not
    always optimumly placed in the limit of small numbers. Takes no account
    of source surveying time, merely creates a grid on the sky, and follows
    that grid.

    Args:
        n_gen (int): Number of pointings.
        **kwargs (dict): Keyword arguments of go.in_region

    Returns:
        ra, dec: Numpy arrays with coordinates of pointings

    """
    def sample(n=1000, random=True):
        """Length of sunflower-like chain to sample."""
        indices = np.arange(0, n, dtype=float) + 0.5

        # Generate coordinates using golden ratio
        ra = (np.pi * (1 + 5**0.5) * indices) % (2*np.pi)
        dec = np.arccos(1 - 2*indices/n) - 0.5*np.pi

        if random:
            # Start from a random point rather than the south pole
            phi = np.random.uniform(-np.pi, np.pi, 1)  # Random ra
            theta = np.random.uniform(-np.pi/2, np.pi/2, 1)  # Random dec

            # Shift ra from -180 to 180
            ra[ra > np.pi] -= 2*np.pi

            # Save on line length
            sin = np.sin
            cos = np.cos
            arcsin = np.arcsin

            y = sin(ra)
            x = np.tan(dec)*sin(theta) + cos(ra)*cos(theta)
            lon = np.arctan2(y, x) - phi
            lat = arcsin(cos(theta)*sin(dec)-cos(ra)*sin(theta)*cos(dec))

            # Shift ra back to 0-360
            ra[ra < 0] += 2*np.pi
            ra = lon % (2*np.pi)
            dec = lat

        # To degrees
        ra = np.rad2deg(ra)
        dec = np.rad2deg(dec)

        return ra, dec

    def accept(ra, dec):
        gl, gb = go.radec_to_lb(ra, dec, frac=True)
        return go.in_region(ra=ra, dec=dec, gl=gl, gb=gb, **kwargs)

    # Accept-Reject
    n = n_gen
    ra, dec = sample(n)
    mask = accept(ra, dec)

    # While there are not enough valid points, keep generating
    while sum(mask) < n_gen:
        n += (n_gen - sum(mask))
        ra, dec = sample(n)
        mask = accept(ra, dec)

    # Only select valid points
    ra = ra[mask]
    dec = dec[mask]

    # Evenly sample arrays
    idx = np.round(np.linspace(0, len(ra) - 1, n_gen)).astype(int)
    ra = ra[idx]
    dec = dec[idx]

    return ra, dec


if __name__ == '__main__':
    """Plot pointings to test pointing generation."""
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from frbpoppy import Survey

    N_POINTS = 1000

    def plot_coordinates(ra, dec):
        """Plot coordinate in 3D plot."""
        dec = np.deg2rad(dec)
        ra = np.deg2rad(ra)

        x = np.cos(ra) * np.cos(dec)
        y = np.sin(ra) * np.cos(dec)
        z = np.sin(dec)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection=Axes3D.name)
        p = ax.scatter(x, y, z, c=np.arange(0, x.size))
        plt.colorbar(p)
        ax.axes.set_xlim3d(left=-1, right=1)
        ax.axes.set_ylim3d(bottom=-1, top=1)
        ax.axes.set_zlim3d(bottom=-1, top=1)
        plt.show()

    _transit = Survey('chime-frb')
    _transit.set_pointings(mount_type='transit', n_pointings=N_POINTS)
    _transit.gen_pointings()
    plot_coordinates(*_transit.pointings)

    _tracking = Survey('perfect-small')
    _tracking.set_pointings(mount_type='tracking', n_pointings=N_POINTS)
    _tracking.gen_pointings()
    plot_coordinates(*_tracking.pointings)
