"""Define class for plotting various parameter distributions."""
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from frbpoppy import CosmicPopulation, Frbcat, unpickle


class PlotParms(object):
    """Plot observed parameters per survey."""

    def __init__(self):
        """Initializing."""
        self.parms = ['fluence', 'dm', 's_peak', 'w_eff']
        self.n_bins = 50
        # If one survey is ever needed, let the default be HTRU
        self.survey = 'HTRU'

    def generate(self):
        """Generate population."""
        self.pop = CosmicPopulation(5000*28,
                                    days=28,
                                    name='standard')

    def read(self):
        """Load population."""
        self.pop = unpickle('standard')

    def preplot(self, parameter):
        """Set plot details.

        Args:
            parameter (type): The parameter to plot.
        """
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.cum = False

        if parameter == 'fluence':
            title = r'Fluence (Jy ms)'
            plt.xscale('log')
            self.logbins = True
        elif parameter == 's_peak':
            title = r'$S_{\text{peak}}$ (Jy)'
            plt.xscale('log')
            plt.yscale('log')
            self.logbins = True
            self.cum = False
        elif parameter == 'dm':
            title = r'Dispersion Measure (\si{\parsec\per\cm\cubed})'
            self.logbins = False
        elif parameter == 'w_eff':
            title = r'Pulse Width (\si{\milli\second})'
            plt.xscale('log')
            plt.yscale('log')
            self.logbins = True
        else:
            self.logbins = False
            title = parameter

        plt.xlabel(title)
        plt.ylabel(r'\# FRB Detections')
        plt.tight_layout()

    def postplot(self, parameter):
        """Post plotting details."""
        # Create new legend handles but use the colors from the existing ones
        handles, labels = self.ax.get_legend_handles_labels()
        new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]
        plt.legend(handles=new_handles, labels=labels)

    def save(self, parameter, plot_type):
        """Save the plot."""
        plt.savefig(f'plots/{parameter}_per_{plot_type}.pdf')

    def gen_bins(self, ps):
        """Generate bins."""
        if self.logbins:
            # Calculate the min and max powers:
            start_power = np.floor(np.log10(min(ps)))
            end_power = np.ceil(np.log10(max(ps)))
            # Generate a range of delimiters in log space
            self.bins = np.logspace(start_power, end_power,
                                    self.n_bins, base=10)
        else:
            self.bins = self.n_bins

    def plot_frbcat(self, parameter):
        """Add a frbcat plot."""
        cat = Frbcat().df
        cat_par = cat[(cat.survey == self.survey)][parameter]

        if len(cat_par) == 0:
            return

        n, bins, patches = self.ax.hist(cat_par,
                                        bins=self.bins,
                                        density=False,
                                        label=f'{self.survey} (frbcat)',
                                        cumulative=self.cum,
                                        histtype='step')
