"""Plot various parameter distributions for different number sidelobes."""
import numpy as np

from frbpoppy import Survey, SurveyPopulation
from parms import PlotParms

MAKE = False


class PlotSidelobe(PlotParms):
    """Plot paramters over sidelobes."""

    def __init__(self):
        """Initializing."""
        PlotParms.__init__(self)
        self.sidelobes = [0, 1, 2, 8]

    def plot(self):
        """Generate plot over surveys."""
        for p in self.parms:

            self.preplot(p)

            for sidelobe in reversed(self.sidelobes):

                survey = Survey(self.survey,
                                gain_pattern='airy',
                                sidelobes=sidelobe,
                                equal_area=self.sidelobes[-1])

                surv_pop = SurveyPopulation(self.pop, survey)

                ps = surv_pop.get(p)
                scaling = np.full((len(ps), 1), surv_pop.rates().f_area)

                self.gen_bins(ps)

                n, bins, patches = self.ax.hist(ps,
                                                bins=self.bins,
                                                density=False,
                                                weights=scaling,
                                                label=f'{sidelobe} sidelobes',
                                                cumulative=self.cum,
                                                histtype='step')

            self.postplot(p)
            self.save(p, 'sidelobe')


if __name__ == '__main__':
    p = PlotSidelobe()

    # Get population
    if MAKE:
        p.generate()
    else:
        p.read()

    p.plot()
