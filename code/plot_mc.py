"""Use Bokeh to plot a population or the results of a Monte Carlo."""

from bokeh.io import output_file, show
from bokeh.layouts import widgetbox
from bokeh.models.widgets import RadioButtonGroup, RangeSlider, Select


class Plot():
    """Class to plot a population or the output of a Monte Carlo."""

    def __init__(self, files=[]):
        """Initializing."""
        self.surveys = ['A', 'B']

        self.pars = {'dec': 'Declination (°)',
                     'dist': 'Distance (Gpc)',
                     'dm_host': 'Dispersion Measure - Host (pc/cm^3)',
                     'dm_igm': 'Dispersion Measure - IGM (pc/cm^3)',
                     'dm_mw': 'Dispersion Measure - Milky Way (pc/cm^3)',
                     'dm': 'Dispersion Measure (pc/cm^3)',
                     'fluence': 'Fluence (Jy*ms)',
                     'gb': 'Galactic Latitude (degrees)',
                     'gl': 'Galactic Longitude (degrees)',
                     'gx': 'Galactic X (Gpc)',
                     'gy': 'Galactic Y (Gpc)',
                     'gz': 'Galactic Z (Gpc)',
                     'lum_bol': 'Luminosity - Bolometric (10^30 ergs/s)',
                     'ra': 'Right Ascension (°)',
                     's_peak': 'Peak Flux Density (Jy)',
                     'snr': 'Signal to Noise Ratio',
                     'w_eff': 'Pulse Width - Effective (ms)',
                     'w_int': 'Pulse Width - Intrinsic (ms)',
                     'z': 'Redshift'}

    def pop(self):
        """Plot a populations."""
        return 'TODO'

    def mc(self):
        """Plot the results of a Monte Carlo run."""
        self.in_pars = ['w_int']
        self.out_pars = ['dec', 'dist', 'dm', 'fluence', 'ra', 'snr', 's_peak',
                         'w_eff', 'z']
        # Output file
        output_file('test.html')

        # Setup survey choice
        sel = ['All'] + self.surveys
        survey_sel = RadioButtonGroup(labels=sel, active=0)

        # Setup observed parameter choice
        opt = [self.pars[o] for o in self.out_pars]
        obs_sel = Select(title="Observed parameter:",
                         value=opt[-2],
                         options=opt)

        # Setup input parameters
        sliders = []
        for p in self.in_pars:
            # TODO convert min and max to values obtained from table
            # r_min = min(table)
            # r_max = max(table)
            # r_step = table[1] - table[0]
            t = self.pars[p]
            r = RangeSlider(start=0, end=10, range=(0,9), step=.1, title=t)
            sliders.append(r)

        show(widgetbox(survey_sel, obs_sel, *sliders))

Plot().mc()
