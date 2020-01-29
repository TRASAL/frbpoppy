"""Plot pointings to test pointing generation."""
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
import numpy as np
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
    ax = fig.add_subplot(111, projection='3d')
    p = ax.scatter(x, y, z, c=np.arange(0, x.size))
    plt.colorbar(p)
    ax.axes.set_xlim3d(left=-1, right=1)
    ax.axes.set_ylim3d(bottom=-1, top=1)
    ax.axes.set_zlim3d(bottom=-1, top=1)
    plt.show()


if __name__ == '__main__':
    transit = Survey('chime')
    transit.set_pointings(mount_type='fixed', n_pointings=N_POINTS)
    transit.gen_pointings()
    plot_coordinates(*transit.pointings)

    tracking = Survey('perfect-small')
    tracking.set_pointings(mount_type='tracking', n_pointings=N_POINTS)
    tracking.gen_pointings()
    plot_coordinates(*tracking.pointings)
