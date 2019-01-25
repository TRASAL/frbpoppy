"""Calculate where null points of an Airy pattern lie."""
from scipy.special import j1
import numpy as np
import matplotlib.pyplot as plt

STEPSIZE = 1e-6
PLOT = True

x_range = np.arange(0, 50, STEPSIZE)
y_range = 4*(j1(x_range)/x_range)**2
nulls = []

# Find where curve changes direction
ind = np.diff(np.sign(np.diff(y_range)))
x_null = x_range[1:-1][ind > 0]
y_null = y_range[1:-1][ind > 0]

print(x_null)

if PLOT:
    title = "Bessel function over x"

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.title(title)
    plt.plot(x_range, y_range)
    plt.scatter(x_null, y_null)
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig('plots/null_sidelobes.pdf')
