"""Define covenience functions."""
import matplotlib.pyplot as plt
import numpy as np
import os

from frbpoppy.log import pprint
import frbpoppy.galacticops as go


def hist(parameter, bin_type='lin', n_bins=25, norm='max', edges=True,
         bins=None):
    """Bin up a parameter either in a lin or log space.

    Why is this not a standard option in numpy or matplotlib?

    Args:
        parameter (array): To be binned
        bin_type (str): Either 'lin' or 'log'
        n_bins (int): Number of bins. Can be overriden internally
        norm (bool): Whether to normalise to 'max' or 'prob' or none

    Returns:
        tuple: bin centers, values per bin

    """
    # Drop NaN-values
    parameter = parameter[~np.isnan(parameter)]

    # Determine number of bins
    if n_bins != 25:
        pass
    elif len(parameter) < 50:
        n_bins = 15
    elif len(parameter) > 500:
        n_bins = 50

    # Determine type of binning
    if bin_type == 'lin':
        _bins = n_bins
    elif bin_type == 'log':
        min_f = np.log10(np.min(parameter[parameter != 0]))
        max_f = np.log10(max(parameter))
        _bins = np.logspace(min_f, max_f, n_bins)

    # Allow for custom bins
    if bins is not None:
        _bins = bins

    # Allow for probability weighting
    weights = None
    if norm == 'prob':
        weights = np.ones(len(parameter)) / len(parameter)

    # Bin
    n, bin_edges = np.histogram(parameter, bins=_bins, weights=weights)

    if norm == 'max':
        n = n/max(n)  # Normalise

    # Centre bins
    bins = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Ensure there are edges on the outer bins of the histograms
    if edges:
        if bin_type == 'lin':
            bin_dif = np.diff(bins)[-1]
            bins = np.insert(bins, 0, bins[0] - bin_dif)
            bins = np.insert(bins, len(bins), bins[-1] + bin_dif)
            n = np.insert(n, 0, 0)
            n = np.insert(n, len(n), 0)
        else:
            bin_dif = np.diff(np.log10(bins))[-1]
            bins = np.insert(bins, 0, 10**(np.log10(bins[0])-bin_dif))
            bins = np.insert(bins, len(bins), 10**(np.log10(bins[-1])+bin_dif))
            n = np.insert(n, 0, 0)
            n = np.insert(n, len(n), 0)

    return bins, n


def plot_aa_style(cols=1):
    """Use a plotting style specifically for Astronomy & Astrophyiscs."""
    # Use a style file
    dir = os.path.abspath(os.path.dirname(__file__))
    style_file = os.path.join(dir, './aa.mplstyle')
    plt.style.use(style_file)

    # Check whether interactive backend is present
    if os.name == 'posix' and "DISPLAY" not in os.environ:
        # On a cluster
        plt.switch_backend('Agg')
        # Use just the basic latex package allowing \text to be used
        plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

    if cols == 2:
        plt.rcParams["figure.figsize"] = (5.75373, 3.556)


def rel_path(path):
    """Return the relative path name for this directory."""
    return os.path.join(os.path.abspath(os.path.dirname(__file__)), path)


def limit_ra_dec(pop, pointings):
    """Doesn't work at boundary coordinates."""
    pprint('Limiting coordinates')

    def sample(n_gen):
        u = np.random.uniform
        ra = u(0, 360, n_gen)
        dec = np.rad2deg(np.arccos(u(-1, 1, n_gen))) - 90
        return ra, dec

    def accept(ra, dec):
        coords = np.full(len(ra), False)
        r = np.sqrt(180/np.pi)
        for p_ra, p_dec in pointings:
            limits = go.separation(ra, dec, p_ra, p_dec) < r
            coords[limits] = True
        return coords

    # Limit population to smaller area
    # Sample RA, dec
    ra, dec = sample(pop.n_gen)
    mask = accept(ra, dec)
    reject, = np.where(~mask)
    while reject.size > 0:
        fill_ra, fill_dec = sample(reject.size)
        mask = accept(fill_ra, fill_dec)
        ra[reject[mask]] = fill_ra[mask]
        dec[reject[mask]] = fill_dec[mask]
        reject = reject[~mask]

    frbs = pop.frbs
    frbs.ra = ra
    frbs.dec = dec

    # Convert to galactic coordinates
    frbs.gl, frbs.gb = go.radec_to_lb(frbs.ra, frbs.dec, frac=True)
    pop.gen_gal_coords()

    return pop


def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''
    from mpl_toolkits.mplot3d import Axes3D

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # # The plot bounding box is a sphere in the sense of the infinity
    # # norm, hence I call half the max range the plot radius.
    # plot_radius = 0.5*max([x_range, y_range, z_range])
    #
    # ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    # ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    # ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

    norm = max(x_range, y_range, z_range)
    x_scale = x_range / norm
    y_scale = y_range / norm
    z_scale = z_range / norm
    scale = np.diag([x_scale, y_scale, z_scale, 1.0])
    scale = scale*(1.0/scale.max())
    scale[3,3] = 1.0

    def short_proj():
      return np.dot(Axes3D.get_proj(ax), scale)

    ax.get_proj = short_proj
