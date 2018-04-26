"""Test whether beam patterns work."""
import numpy as np
from bokeh.plotting import figure, show
from bokeh.palettes import Category10
from bokeh.models import Range1d, Span
import itertools

from frbpoppy.survey import Survey

# 'gaussian', 'airy', 'tophat', 'perfect'
PATTERNS = ['gaussian', 'airy', 'parkes']
# 'beam_pattern' or 'hist_int_pro'
PLOT = 'hist_int_pro'

n = 100000

# Set up plot
if PLOT == 'beam_pattern':
    p = figure()
elif PLOT == 'hist_int_pro':
    p = figure(x_axis_type='log')


# Get a colour generator
def color_gen():
    yield from itertools.cycle(Category10[10])
color = color_gen()

for pattern in PATTERNS:
    s = Survey('APERTIF', gain_pattern=pattern)

    # Plot beam pattern
    if PLOT == 'beam_pattern':
        data = [s.intensity_profile(test=True, sidelobe=1) for e in range(n)]
        x, y = zip(*data)
        print(x, y)
        p.circle(x, y, legend=pattern, color=next(color))

    elif PLOT == 'hist_int_pro':
        data = [s.intensity_profile(test=True, sidelobe=1) for e in range(n)]
        _, ints = zip(*data)

        ints = [i for i in ints if i != 0]

        # Bin up
        log_bins = np.logspace(np.log10(min(ints)), np.log10(max(ints)), 500)
        hist, edges = np.histogram(ints, bins=log_bins)
        norm = hist / max(hist)

        p.xaxis.axis_label = 'Intensity Profile Factor'
        p.yaxis.axis_label = 'Fraction'

        # Find optimum lower limit
        bottom = 10**(np.log10(1/n) - 1)

        p.quad(top=norm,
               bottom=bottom,
               left=edges[:-1],
               right=edges[1:],
               alpha=0.5,
               legend=pattern,
               color=next(color))


# Add beam size
if PLOT == 'beam_pattern':
    radius = s.fwhm / 2

    # Plot vertical line
    vline = Span(location=radius,
                 dimension='height',
                 line_color='red',
                 line_width=2)

    p.renderers.extend([vline])

show(p)
