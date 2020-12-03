"""Scripts to understand what inputs mean for a lognormal distribution."""
import numpy as np
import matplotlib.pyplot as plt
from frbpoppy import hist


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


fig, axes = plt.subplots(1, 1, sharey=True)

mean = 2.71**2
std = 5
size = int(1e4)

axes.set_xscale('log', basex=2.71)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

data = lognormal(mean, std, size)
xy = hist(data, bin_type='log')
label = 'Distribution with input (frbpoppy)'
axes.step(*xy, where='mid', color=colors[0], label=label)
axes.axvline(mean, color=colors[0], linestyle='dotted')
axes.text(mean, 1, mean)
axes.axvline(mean+std, color=colors[0], linestyle='dotted')
axes.text(mean+std, 1, f'{mean}+{std}')


data = np.random.lognormal(mean, std, size)
xy = hist(data, bin_type='log')
label = 'Underlying normal input (numpy)'
axes.step(*xy, where='mid', color=colors[1], label=label)
axes.axvline(mean, color=colors[1], linestyle='dotted')
axes.text(mean, 1, mean)
axes.axvline(mean+std, color=colors[1], linestyle='dotted')
axes.text(mean+std, 1, f'{mean}+{std}')

plt.legend()
plt.show()
