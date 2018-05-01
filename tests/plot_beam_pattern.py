"""Test whether beam patterns work."""
import numpy as np
from bokeh.plotting import figure, show
from bokeh.palettes import Category10
from bokeh.models import Span
import itertools

from frbpoppy.survey import Survey

PATTERNS = ['airy']  # 'gaussian', 'airy', 'tophat', 'perfect'
PLOT = 'hist_int_pro'  # 'beam_pattern' or 'hist_int_pro'
EQUAL_AREA = True
DIMENSIONS = 1

n = 100000

# Set up plot
if PLOT == 'beam_pattern':
    p = figure(y_axis_type='log')
elif PLOT == 'hist_int_pro':
    p = figure(x_axis_type='log')

# Get a sequence of colours
def color_gen():
    """Generate a sequence of colours."""
    yield from itertools.cycle(Category10[10])
color = color_gen()

# Generate a data set
def get_int_pro(sidelobes=1):
    """Get a whole range of intensity profiles."""
    args = {'test': True,
            'sidelobes': sidelobes,
            'equal_area': EQUAL_AREA,
            'dimensions': DIMENSIONS}
    data = [s.intensity_profile(**args) for e in range(n)]
    offset, int_pro = zip(*data)
    return offset, int_pro


def plot_hist(ints, pattern='', i=None):
    """Plot a histogram of a list."""
    ints = [i for i in ints if i != 0]

    # Bin up
    log_bins = np.logspace(np.log10(min(ints)), np.log10(max(ints)), 500)
    hist, edges = np.histogram(ints, bins=log_bins)
    norm = hist / max(hist)

    p.xaxis.axis_label = 'Intensity Profile Factor'
    p.yaxis.axis_label = 'Fraction'

    # Find optimum lower limit
    bottom = 10**(np.log10(1/n) - 1)

    if i is not None:
        pattern += f'_{i}'

    p.quad(top=hist,
           bottom=bottom,
           left=edges[:-1],
           right=edges[1:],
           alpha=0.5,
           legend=pattern,
           color=next(color))


for pattern in PATTERNS:
    s = Survey('APERTIF', gain_pattern=pattern)

    # Plot beam pattern
    if PLOT == 'beam_pattern':

        # If wanting to compare sidelobes
        if len(PATTERNS) == 1 and pattern == 'airy':
            for i in range(3, -1, -1):
                x, y = get_int_pro(sidelobes=i)
                p.circle(x, y, legend=f'{pattern}_{i}', color=next(color))

        else:
            x, y = get_int_pro()
            p.circle(x, y, legend=pattern, color=next(color))

    elif PLOT == 'hist_int_pro':

        # If wanting to compare sidelobes
        if len(PATTERNS) == 1 and pattern == 'airy':
            for i in range(3, -1, -1):
                _, ints = get_int_pro(sidelobes=i)
                plot_hist(ints, pattern=pattern, i=i)
        else:
            _, ints = get_int_pro()
            plot_hist(ints, pattern=pattern)


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
