import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from copy import deepcopy

from frbpoppy import CosmicPopulation, SurveyPopulation, Survey
from frbpoppy import poisson_interval
from tests.convenience import plot_aa_style, rel_path

ACTUAL_RATES = False
SURVEYS = ['askap-fly', 'fast-crafts', 'chime-frb', 'wrst-apertif']
SIZE = 1e5


class Rates:

    def __init__(self, size, survey_names, rates=None):
        self.size = size
        self.survey_names = survey_names
        self.surveys = [self.gen_survey(s) for s in survey_names]
        self.data = []
        self.rates = rates

        self.zs = [0.01, 2, 6]

        # Generate populations
        self.gen_def_pop()
        self.gen_cos_pops()
        self.gen_lum_pops()
        self.gen_si_pops()
        self.gen_w_pops()
        self.df = pd.DataFrame(self.data)
        self.plot()

    def loop_surveys(self, cosmic_pop, z_max):
        for survey in self.surveys:
            surv_pop = SurveyPopulation(cosmic_pop, survey)
            surv_pop.z_max = z_max
            rate = surv_pop.source_rate.det
            n = surv_pop.n_sources()
            errs = [np.abs(e-n) for e in poisson_interval(n)]
            d = {'survey': survey.name,
                 'pop': cosmic_pop.name,
                 'rate': rate,
                 'err_low': errs[0],
                 'err_high': errs[1],
                 'z_max': z_max}
            self.data.append(d)

    def gen_def_pop(self):
        self.default_pop = CosmicPopulation.simple(self.size)
        self.default_pop.set_lum(model='constant', value=1e43)
        self.default_pop.generate()

    def gen_survey(self, name):
        survey = Survey(name)
        survey.set_beam('gaussian')
        return survey

    def gen_cos_pops(self):
        cosmo_pop = deepcopy(self.default_pop)
        for z in self.zs:
            cosmo_pop.set_dist(z_max=z)
            cosmo_pop.name = r'z$_{\text{max}}$=' + str(z)
            cosmo_pop.generate()
            self.loop_surveys(cosmo_pop, z)

    def gen_lum_pops(self):
        lum_pop = deepcopy(self.default_pop)
        for z in [0.01, 6]:
            # Standard candles
            lum_pop.set_lum(model='constant', value=1e40)
            lum_pop.generate()
            lum_pop.name = f'std candle'
            self.loop_surveys(lum_pop, z)

            # Powerlaw with slope
            power = -1
            lum_pop.set_lum(model='powerlaw', low=1e40, high=1e43, power=power)
            lum_pop.generate()
            lum_pop.name = f'li={power}'
            self.loop_surveys(lum_pop, z)

            # Powerlaw with slope
            power = -2
            lum_pop.set_lum(model='powerlaw', low=1e40, high=1e43, power=power)
            lum_pop.generate()
            lum_pop.name = f'li={power}'
            self.loop_surveys(lum_pop, z)

    def gen_si_pops(self):
        si_pop = deepcopy(self.default_pop)
        for z in [0.01, 6]:
            si_pop.set_dist(z_max=z)
            for si in [-2, 0, 2]:
                si_pop.set_si(model='constant', value=si)
                si_pop.name = f'si={si}'
                si_pop.generate()
                self.loop_surveys(si_pop, z)

    def gen_w_pops(self):
        w_pop = deepcopy(self.default_pop)

        for z in [0.01, 6]:
            w_pop.set_dist(z_max=z)

            # Constant
            w_pop.set_w(model='constant', value=10)
            w_pop.generate()
            w_pop.name = f'constant'
            self.loop_surveys(w_pop, z)

            # Normal
            w_pop.set_w(model='gauss', mean=10, std=10)
            w_pop.generate()
            w_pop.name = f'normal'
            self.loop_surveys(w_pop, z)

            # Lognormal
            w_pop.set_w(model='lognormal', mean=10, std=10)
            w_pop.generate()
            w_pop.name = f'lognormal'
            self.loop_surveys(w_pop, z)

    def plot(self):
        plot_aa_style()
        plt.rcParams["figure.figsize"] = (5.75373, 5.75373)
        f, self.axes = plt.subplots(2, 2, sharex='col', sharey='row')

        self.linestyles = ['solid', 'dashed', 'dotted']
        self.colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
        # Matching redshifts to linestyles
        self.zs = self.df.z_max.unique()
        self.lz = dict(zip(self.zs, self.linestyles))

        # Set vertical spacing
        self.internal_survey_spacing = 0.1
        self.external_survey_spacing = 0.2

        for ax in self.axes.flat:
            ax.set_xscale('log')
            ax.set_xlim(1e-6, 1e6)

        self.plot_cosmo()
        self.plot_lum()
        self.plot_si()
        self.plot_w()

        plt.tight_layout()
        p = f'plots/flattening_rates.pdf'
        plt.savefig(rel_path(p))

    def plot_rate(self, ax, df, subplot=None):
        y = 0
        label_ys = []
        surveys = df.survey.unique()
        norm_rate = df.iloc[0].rate

        if subplot == 'cosmo':
            y -= 1.5*self.internal_survey_spacing

        for survey, group in df.groupby('survey', sort=False):
            y_min = y
            color_ix = 0
            for pop, group in group.groupby('pop', sort=False):

                # Set colour per population type
                colour = self.colours[color_ix]
                if subplot == 'cosmo':
                    colour = 'grey'
                color_ix += 1
                ls_ix = 0
                for _, row in group.iterrows():

                    # Set linestyle per population
                    ls = self.lz[row.z_max]
                    ls_ix += 1
                    label = pop

                    # Set up legend for just the first set
                    if subplot == 'cosmo':
                        if survey != surveys[0]:
                            label = f'_{label}'
                    else:
                        if ls_ix > 1 or survey != surveys[0]:
                            label = f'_{label}'

                    # Plot errorbars
                    errs = [row.err_low, row.err_high]
                    line = ax.errorbar(row.rate/norm_rate, y,
                                       xerr=np.array([errs]).T,
                                       fmt='x',
                                       color=colour,
                                       label=rf'{label}')
                    line[-1][0].set_linestyle(ls)

                    # Shift vertical position to next line
                    y -= self.internal_survey_spacing

            # Show survey interval
            # overflow = 0.25*self.external_survey_spacing
            # y_max = y + self.internal_survey_spacing
            # ax.plot([ax.get_xlim()[0]]*2,
            #         [y_min+overflow, y_max-overflow],
            #         color='k', lw=4)

            if ACTUAL_RATES:  # Add line with actual rate
                ax.plot([self.rates[survey]]*2,
                        [y_min, y + self.internal_survey_spacing],
                        color='grey')

            # Calculate coordinates for y axis
            y_width = [y_min, y+self.internal_survey_spacing]
            label_ys.append(np.mean(y_width))
            y -= self.external_survey_spacing

            if subplot == 'cosmo':
                y -= 3*self.internal_survey_spacing

        ax.set_yticks(label_ys)
        ax.set_yticklabels(surveys)
        ax.legend(prop={'size': 8}, markerscale=0)

    def plot_cosmo(self):
        ax = self.axes[0, 0]
        ax.set_title('Cosmology')
        filtered_df = self.df[(self.df['pop'].str.startswith('z$_'))]
        self.plot_rate(ax, filtered_df, subplot='cosmo')
        ax.invert_yaxis()

    def plot_lum(self):
        ax = self.axes[0, 1]
        ax.set_title('Luminosity')
        p = self.df['pop']
        filtered_df = self.df[((p == 'std candle') | (p.str.startswith('li')))]
        self.plot_rate(ax, filtered_df)
        ax.invert_yaxis()

    def plot_si(self):
        ax = self.axes[1, 0]
        ax.set_xlabel(r'Rate (day$^{-1}$)')
        ax.set_title('Spectral index')
        self.plot_rate(ax, self.df[(self.df['pop'].str.startswith('si'))])
        ax.invert_yaxis()

    def plot_w(self):
        ax = self.axes[1, 1]
        ax.set_xlabel(r'Rate (day$^{-1}$)')
        ax.set_title('Pulse width')
        keep_list = ['constant', 'lognormal', 'normal']
        self.plot_rate(ax, self.df[self.df['pop'].isin(keep_list)])
        ax.invert_yaxis()


if __name__ == '__main__':
    rates = {'wrst-apertif': 1, 'parkes-htru': 2, 'fast-crafts': 5, 'askap': 10}
    Rates(SIZE, SURVEYS, rates)
