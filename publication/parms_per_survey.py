"""Plot various parameter distribution for different surveys."""
import numpy as np

from frbpoppy import Survey, SurveyPopulation
from parms import PlotParms

MAKE = False
FRBCAT = True


class PlotSurveys(PlotParms):
    """Plot paramters over surveys."""

    def __init__(self):
        """Initializing."""
        PlotParms.__init__(self)
        self.patterns = ['parkes', 'apertif']
        self.surveys = ['HTRU', 'APERTIF']

    def plot(self):
        """Generate plot over surveys."""
        for p in self.parms:

            self.preplot(p)

            for i, pattern in enumerate(self.patterns):
                self.survey = self.surveys[i]
                survey = Survey(self.survey, gain_pattern=pattern)
                surv_pop = SurveyPopulation(self.pop, survey)
                ps = surv_pop.get(p)
                scaling = np.full((len(ps), 1), surv_pop.rates().f_area)

                self.gen_bins(ps)

                n, bins, patches = self.ax.hist(ps,
                                                bins=self.bins,
                                                density=False,
                                                weights=scaling,
                                                label=self.survey,
                                                cumulative=self.cum,
                                                histtype='step')

                if FRBCAT:
                    self.plot_frbcat(p)

            self.postplot(p)
            self.save(p, 'survey')


if __name__ == '__main__':
    p = PlotSurveys()
    # Get population
    if MAKE:
        p.generate()
    else:
        p.read()

    p.plot()
