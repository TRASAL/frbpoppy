"""Plot how rates vary with luminosity index."""
from collections import namedtuple
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, paths

MAKE = False


class Rates:
    """Hold various functions to be able to plot rates per luminosity index."""

    def __init__(self):
        """Initializing."""
        self.surveys = ('htru', 'apertif', 'utmost', 'askap-fly')
        self.lis = np.linspace(-3, 1, num=5)  # Luminosity Indices

        # Data Point
        self.DP = namedtuple('Data', ['li', 'survey', 'exp', 'fluence', 'det'])
        self.data = []
        self.path = f'{paths.results()}rates_li.csv'

    def generate(self):
        """Generate data."""
        for li in self.lis:
            n = 5000
            days = 28
            pop = CosmicPopulation(n*days, days=days, lum_index=li)

            for survey in self.surveys:
                surv = Survey(survey, gain_pattern='airy')
                surv_pop = SurveyPopulation(pop, surv)

                exp = surv_pop.rates().exp
                det = surv_pop.rates().det
                try:
                    fluence = min(surv_pop.get('fluence'))
                except ValueError:
                    fluence = float('Nan')

                self.data.append(self.DP(li, survey, exp, fluence, det))

        self.df = pd.DataFrame.from_records(self.data, columns=self.DP._fields)
        self.save()

    def save(self):
        """Save to CSV."""
        df = pd.DataFrame.from_records(self.data, columns=self.DP._fields)
        df.to_csv(self.path)

    def read(self):
        """Import data from CSV."""
        self.data = []
        self.df = pd.read_csv(self.path)

        for i, r in self.df.iterrows():
            self.data.append(self.DP(r.li, r.survey, r.exp, r.fluence, r.det))

    def plot(self):
        """Plot the expected rate."""
        # Create plot
        fig, ax = plt.subplots()

        # Do stuff with data
        self.df.sort_values('li', inplace=True)
        for survey, group in self.df.groupby('survey'):
            plt.plot(group.li, group.exp, 'o-', label=survey)

        ax.set_xlabel('Luminosity Index')
        ax.set_ylabel('Days per FRB')
        plt.legend()
        plt.tight_layout()
        plt.yscale('log')
        plt.savefig('plots/rates_per_lum_index.pdf')

        # Create plot with fluence limit
        fig, ax = plt.subplots()

        # Do stuff with data
        self.df.sort_values('li', inplace=True)
        for survey, group in self.df.groupby('survey'):
            plt.plot(group.li, group.fluence, 'o-', label=survey)

        ax.set_xlabel('Luminosity Index')
        ax.set_ylabel('Fluence Limit (Jy ms)')
        plt.legend()
        plt.tight_layout()
        plt.yscale('log')
        plt.savefig('plots/rates_per_lum_index_fluence_limit.pdf')

if __name__ == '__main__':
    r = Rates()
    if MAKE:
        r.generate()
    else:
        r.read()
    r.plot()
