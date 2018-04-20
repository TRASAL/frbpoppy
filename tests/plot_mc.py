"""
Plot the results of a Monte Carlo simulation.

Run with::

    $ bokeh serve --show tests/plot_mc.py

"""
from frbpoppy.plot_mc import Plot

Plot().mc()

if __name__ == '__main__':
    print('Please launch with: bokeh serve --show tests/plot_mc.py')
