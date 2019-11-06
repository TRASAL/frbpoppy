"""Test the behaviour of clustered burts."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.special import gamma
import scipy.special as special
from scipy.optimize import curve_fit
from scipy.stats import norm

from convenience import plot_aa_style, rel_path

M = 700
DAYS = 50
N_FRBS = int(1e5)

days = np.arange(1, DAYS, 1)
n_bursts = np.arange(0, M, 1)


def gen_clustered_times(n, r=5.7, k=0.34, max_days=0):
    """Generate burst times following Oppermann & Pen (2017)."""
    # Determine the maximum possible number of bursts per source to include
    # m = int(round(np.log10(self.n_gen))*2) - 1
    m = M
    if m < 1:
        m = 1
    dims = (int(n), m)

    lam = 1/(r*gamma(1 + 1/k))
    time = lam*np.random.weibull(k, dims).astype(np.float32)
    time = np.cumsum(time, axis=1)   # This is in fraction of days

    # The first burst time is actually the time since the previous one
    # You want to be at in a random time in between those
    time_offset = np.random.uniform(0, time[:, 0])
    time -= time_offset[:, np.newaxis]

    return time


def get_probabilities(normalise=False):
    """Get the number of bursts per maximum time.

    Args:
        normalise (type): Whether to normalise at each maximum time, such that
            there's a number of bursts at each time stamp having a maximum
            probability of one.

    Returns:
        array: Number of bursts per maximum time.

    """
    n_frbs = N_FRBS
    days = np.arange(1, DAYS, 1)
    n_bursts = np.arange(0, M, 1)
    xx, yy = np.meshgrid(n_bursts, days)
    prob = np.full((len(days), len(n_bursts)), np.nan)

    # Mask any frbs over the maximum time
    time = gen_clustered_times(n_frbs)

    for i, day in enumerate(days):
        t = np.copy(time)
        t[(t > day)] = np.nan
        m_bursts = (~np.isnan(t)).sum(1)
        unique, counts = np.unique(m_bursts, return_counts=True)
        prob[(i, unique)] = counts

    prob = prob.T/n_frbs

    if normalise:
        # Normalise at each maximum time (highest chance becomes one)
        prob = prob / np.nanmax(prob, axis=0)

    return prob


def plot(prob):
    """Plot the number of bursts seen over a maximum time."""
    plot_aa_style(cols=2)

    prob[np.isnan(prob)] = 0

    fig, ax = plt.subplots()
    extent = [min(days)-0.5, max(days)+0.5, min(n_bursts), max(n_bursts)]
    plt.imshow(prob, norm=LogNorm(), origin='lower', aspect='auto',
               interpolation='none', extent=extent)
    plt.colorbar(orientation='vertical')
    ax.set_title('Probability of M bursts within time')
    ax.set_xlabel('Maximum time (days)')
    ax.set_ylabel('M bursts')
    plt.tight_layout()
    plt.savefig(rel_path(f'./plots/prob_m_bursts.pdf'))


def calc_m(max_time, n_gen):
    """Calculate the number of expected bursts (Oppermann clustering).

    Note this function only works if following the exact shape parameters of
    the Oppermann and Pen paper.

    Args:
        max_time (float): Maximum amount of time spent on .
        n_gen (type): Description of parameter `n_gen`.

    Returns:
        type: Description of returned object.

    """
    mu = 5.8*max_time + 4.86
    sigma = -0.014*max_time**2+1.80*max_time+9.34
    prob = 1/(sigma*np.sqrt(2*np.pi)*n_gen)
    discrim = -2 * sigma**2 * np.log(prob * sigma * np.sqrt(2*np.pi))
    return mu + np.sqrt(discrim)


def plot_fitted_functions_too(prob):
    """Plot the fitted function over the probabilities."""
    prob[np.isnan(prob)] = 0

    for n_gen in (1, 1e1, 1e2, 1e3, 1e4, 1e5):
        xs = np.arange(1, prob.shape[1]+1)
        ys = calc_m(xs, n_gen)
        plt.plot(xs, ys, label=f'n = 1e{int(np.log10(n_gen))}')

    plt.legend()
    plt.tight_layout()

    plt.savefig(rel_path(f'./plots/prob_w_functions.pdf'))


def fit_prob(prob, plot=False):
    """Fit the generated probabilities.

    Honestly, I would prefer it somebody could derive this analytically.
    """
    # prob[np.isnan(prob)] = 0

    # First fit a line along the maximums
    def linear(x, a, b):
        return a*x + b

    maxs = np.nanargmax(prob, axis=0)
    a, b = curve_fit(linear, days, maxs)[0]

    # Then fit gaussians along each point of the line
    def to_fit_gaussian(mn, std):
        m = mn[:, 0]
        n = mn[:, 1][0]
        return norm(a*n + b, std).pdf(m)

    stds = []
    ns = np.arange(1, prob.shape[1]+1)
    for n in ns:
        y = prob[:, n-1]
        x = np.arange(len(y))

        mn = np.c_[x, np.ones(len(x))*n]
        std = curve_fit(to_fit_gaussian, mn, y)[0][0]
        stds.append(std)

        # Plot the histogram.
        if plot:
            plt.plot(x, y, alpha=0.6, color='g')

            # Plot the PDF.
            xmin, xmax = plt.xlim()
            x = np.linspace(xmin, xmax, 100)
            plt.plot(x, norm(a*n + b, std).pdf(x), 'k', linewidth=2)
            title = "Fit results: std = %.2f" % (std)
            plt.title(title)
            plt.show()
            plt.clf()

    # Then fit all of the standard deviations with a quadratic formula
    def to_fit_std(n, c, d, e):
        return c*n**2+d*n+e

    c, d, e = curve_fit(to_fit_std, ns, stds)[0]

    return a, b, c, d, e


if __name__ == '__main__':
    prob = get_probabilities(normalise=True)
    plot(prob)
    a, b, c, d, e = fit_prob(prob)
    m = 'The number of bursts can be described by: \n'
    m += f'mu = {a:.2f}*max_time + {b:.2f}\n'
    m += f'sigma = {c:.2f}*max_time**2+{d:.2f}*max_time+{e:.2f} \n'
    m += 'prob = 1/(sigma*np.sqrt(2*np.pi)*n_gen) \n'
    m += 'discrim = -2 * sigma**2 * np.log(prob * sigma * np.sqrt(2*np.pi)) \n'
    m += 'm = mu + np.sqrt(discrim) \n'
    print(m)
