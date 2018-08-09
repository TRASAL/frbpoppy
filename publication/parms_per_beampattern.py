"""Plot various parameter distributions for different beam types."""
import numpy as np

from frbpoppy import Survey, SurveyPopulation
from parms import PlotParms

MAKE = False


class PlotBeamPattern(PlotParms):
    """Plot paramters over beam patterns."""

    def __init__(self):
        """Initializing."""
        PlotParms.__init__(self)
        self.patterns = ['perfect', 'tophat', 'gaussian', 'airy']

    def plot(self):
        """Generate plot over surveys."""
        for p in self.parms:

            self.preplot(p)

            for pattern in self.patterns:

                survey = Survey(self.survey, gain_pattern=pattern)
                surv_pop = SurveyPopulation(self.pop, survey)
                ps = surv_pop.get(p)
                scaling = np.full((len(ps), 1), surv_pop.rates().f_area)

                self.gen_bins(ps)

                n, bins, patches = self.ax.hist(ps,
                                                bins=self.bins,
                                                density=False,
                                                weights=scaling,
                                                label=pattern,
                                                cumulative=self.cum,
                                                histtype='step')

            self.postplot(p)
            self.save(p, 'beampattern')


if __name__ == '__main__':
    p = PlotBeamPattern()

    # Get population
    if MAKE:
        p.generate()
    else:
        p.read()

    p.plot()
