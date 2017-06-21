"""Use Bokeh to plot a population or the results of a Monte Carlo."""

from bokeh.io import output_file, show
from bokeh.layouts import widgetbox
from bokeh.models.widgets import RadioButtonGroup, RangeSlider, Select


class Plot():
    """Class to plot a population or the output of a Monte Carlo."""

    def __init__(self, files=[]):
        """Initializing."""
        self.surveys = ['A', 'B']
        self.obs_par = ['eff_wid', 'obs_test']

    def pop(self):
        """Plot a populations."""
        return 'TODO'

    def mc(self):
        """Plot the results of a Monte Carlo run."""

        output_file("test.html")

        # Setup survey choice
        sel = ['All'] + self.surveys
        survey_sel = RadioButtonGroup(labels=sel, active=0)

        # Setup observed parameter choice
        obs_sel = Select(title="Observed parameter:", options=self.obs_par)

        # Setup input parameters
        # TODO turn into dictionary loop
        t = 'Pulse Width - Intrinsic (ms)'
        w_min_sel = RangeSlider(start=0, end=10, range=(0,9), step=.1, title=t)


        show(widgetbox(survey_sel, obs_sel, w_min_sel))

Plot().mc()
