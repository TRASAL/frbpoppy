"""Calculate where null points of an Airy pattern lie."""
from scipy.special import j1
from numpy import arange
from bokeh.plotting import figure, show

STEPSIZE = 1e-6
PLOT = False

x_range = arange(0, 100, STEPSIZE)
y_range = 4*(j1(x_range)/x_range)**2
nulls = []

for i in range(2, len(x_range)):
    if y_range[i-2] > y_range[i-1] < y_range[i]:
        null = (x_range[i-1], y_range[i-1])
        nulls.append(null)

x_null, y_null = zip(*nulls)

print(x_null)

if PLOT:
    p = figure(title="Bessel function over x")
    p.line(x_range, y_range)
    p.circle(x_null, y_null)
    show(p)
