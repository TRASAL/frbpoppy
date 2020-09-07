"""Run a Monte Carlo determining best fit parameters.

TODO: Work in progress."""
from frbpoppy import CosmicPopulation, Survey, SurveyPopulation
from frbpoppy import pprint, unpickle, Frbcat, poisson_interval
from joblib import Parallel, delayed
from matplotlib.colors import LogNorm
from scipy.stats import ks_2samp
from tqdm import tqdm
import frbpoppy.paths
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import os
import pandas as pd

from tests.rates.alpha_real import EXPECTED
from tests.convenience import plot_aa_style, rel_path

SURVEY_NAMES = ['apertif', 'askap-fly', 'htru', 'chime']
ALPHAS = np.linspace(-2, -0.5, 3)
LIS = np.linspace(-2, 0, 5)
SIS = np.linspace(-2, 2, 5)
POP_SIZE = 1e6
MAKE = False
TEST_DATA = False
PLOT = True


class Run:
    """Hold information relevant to a single Monte Carlo run."""

    def __init__(self, alpha, li, si):
        self.li = li
        self.si = si
        self.alpha = alpha
        self.survey_name = None
        self.ks_dm = np.nan
        self.ks_snr = np.nan
        self.ks_rate = np.nan
        self.pop_name = None
        self.pop_size = None

    def to_dict(self):
        """Convert to dictionary."""
        return {'survey_name': self.survey_name,
                'ks_dm': self.ks_dm,
                'ks_snr': self.ks_snr,
                'ks_rate': self.ks_rate,
                'li': self.li,
                'si': self.si,
                'alpha': self.alpha,
                'pop_path': self.pop_name}


class MonteCarlo:

    def __init__(self, alphas, lis, sis, survey_names, pop_size):
        self.alphas = alphas
        self.lis = lis
        self.sis = sis
        self.pop_size = pop_size
        self.survey_names = survey_names
        self.set_up_surveys()
        self.set_up_dir()
        self.get_frbcat()
        self.runs = []

    def get_frbcat(self):
        # Only get one-offs
        self.frbcat = Frbcat(repeaters=False, mute=True).df

    def set_up_surveys(self):
        """Set up surveys."""
        self.surveys = []
        for name in self.survey_names:
            survey = Survey(name=name)
            survey.set_beam(model='airy', n_sidelobes=1)
            self.surveys.append(survey)

    def set_up_dir(self):
        """Create subdirectory for saving populations."""
        f = f'{frbpoppy.paths.populations()}mc/'
        if not os.path.isdir(f):
            os.mkdir(f)

    def gen_runs(self, parallel=True):

        def iter_alpha(i, parallel=None):
            runs = []

            alpha = self.alphas[i]
            pop = CosmicPopulation.complex(self.pop_size)
            pop.set_dist(model='vol_co', z_max=1.0, alpha=alpha)
            pop.set_lum(model='constant', value=1)
            pop.generate()

            for li in self.lis:
                pop.set_lum(model='powerlaw', low=1e40, high=1e45, power=li)
                pop.gen_lum()

                for si in self.sis:
                    pop.set_si(model='constant', value=si)
                    pop.gen_si()
                    pop.name = f'mc/complex_alpha_{alpha}_lum_{li}_si_{si}'

                    for survey in self.surveys:
                        surv_pop = SurveyPopulation(pop, survey)
                        surv_pop.save()
                        run = self.pop_to_run(surv_pop)
                        runs.append(run)

            return runs

        if parallel:
            n_cpu = min([3, os.cpu_count() - 1])
            pprint(f'{os.cpu_count()} CPUs available')
            r = range(len(self.alphas))
            p = Parallel(n_jobs=n_cpu)(delayed(iter_alpha)(i) for i in tqdm(r, desc='Parallel processes: '))

        else:
            p = []
            for i in tqdm(range(len(self.alphas)), desc='Alphas'):
                p.append(iter_alpha(i))

        self.runs = [i for sublist in p for i in sublist]

    def pop_to_run(self, surv_pop):
        """Given population information, construct run information."""
        # Save output
        n = surv_pop.name.split('_')
        alpha, li, si = [float(n[i]) for i in (-6, -4, -2)]
        run = Run(alpha, li, si)
        run.survey_name = n[-1]
        run.pop_name = surv_pop.name

        # Add rate details
        sr = surv_pop.source_rate
        run.n_srcs = sr.det
        run.n_days = sr.days
        run.rate = sr.det / sr.days
        rate_errs = [p/sr.days for p in poisson_interval(sr.det, sigma=1)]

        if sr.det == 0:
            return run

        if run.survey_name in EXPECTED:
            exp = EXPECTED[run.survey_name]
        else:
            exp = [np.nan, np.nan]

        actual_rate = exp[0]/exp[1]
        # TODO: Update to include information on the Poisson intervals
        run.ks_rate = np.abs(run.rate - actual_rate)

        mask = (self.frbcat.survey == run.survey_name)
        run.ks_dm = ks_2samp(surv_pop.frbs.dm, self.frbcat[mask].dm)[1]
        run.ks_snr = ks_2samp(surv_pop.frbs.snr, self.frbcat[mask].snr)[1]

        return run

    def get_runs(self):
        """Load in populations to get information on the runs."""
        for alpha in self.alphas:
            for li in self.lis:
                for si in self.sis:
                    for surv in self.survey_names:
                        f = f'mc/complex_alpha_{alpha}_lum_{li}_si_{si}_{surv}'
                        surv_pop = unpickle(f)
                        run = self.pop_to_run(surv_pop)
                        self.runs.append(run)

    def plot(self):
        # Convert to a pandas dataframe
        pprint('Plotting')
        df = pd.DataFrame([r.to_dict() for r in self.runs])

        for ks_type in ('ks_dm', 'ks_snr', 'ks_rate'):
            for i, group in df.groupby('survey_name'):
                self.plot_run(ks_type, group)
            # A combined plot per ks type
            self.plot_run(ks_type, df)
        # nd one total combined plot
        ks_cols = [col for col in df if col.startswith('ks_')]
        df['ks_all'] = df[ks_cols].product(axis=1, skipna=True, min_count=1)
        self.plot_run('ks_all', df)

    def plot_run(self, ks_type, df):
        if TEST_DATA:
            df = self.apply_test(df)

        plot_aa_style()
        plt.rcParams["figure.figsize"] = (5.75373, 5.75373)
        fig, axes = plt.subplots(3, 3, sharex='col')

        if len(df.survey_name.unique()) == 1:
            survey_name = df.survey_name.unique()[0]
        else:
            survey_name = 'combined'

        # Add color grids
        args = {'cmap': 'viridis', 'norm': LogNorm(),
                'vmin': 1e-6, 'vmax': 1e0}
        axes[2, 0].pcolormesh(*self.make_mesh('alpha', 'si', df, ks_type),
                              **args)
        axes[1, 0].pcolormesh(*self.make_mesh('alpha', 'li', df, ks_type),
                              **args)
        im = axes[2, 1].pcolormesh(*self.make_mesh('li', 'si', df, ks_type),
                                   **args)

        # Add histograms
        axes[0, 0].step(*self.make_hist('alpha', df, ks_type), where='mid')
        axes[1, 1].step(*self.make_hist('li', df, ks_type), where='mid')
        axes[2, 2].step(*self.make_hist('si', df, ks_type), where='mid')

        # Add plot information
        axes[0, 2].text(.5, .7, fr'{survey_name} | {ks_type[3:]}',
                        horizontalalignment='center',
                        verticalalignment='center',
                        transform=axes[0, 2].transAxes)

        # Center ticks on pixels
        for i, ax in enumerate(axes.flatten()):
            ax.label_outer()
            x_par = ax.get_xlabel()
            y_par = ax.get_ylabel()
            if x_par:
                xs = np.sort(df[x_par].unique())
                dx = xs[1] - xs[0]
                x2 = xs + dx
                ax.set_xlim(np.min(x2), np.max(x2))
                ax.set_xticks(xs)
            if y_par:
                ys = np.sort(df[y_par].unique())
                y_min = np.min(ys)
                y_max = np.max(ys)
                dy = ys[1] - ys[0]
                y2 = np.linspace(y_min, y_max+dy, dy)-dy/2.
                ax.set_ylim(np.min(y2), np.max(y2))
                ax.set_yticks(y2)

            # Apply plot limits
            for n, p in enumerate(('alpha', 'li', 'si')):
                if i % 3 == n:
                    xs = np.sort(df[p].unique())
                    dx = xs[1] - xs[0]
                    ax.set_xlim(xs[0]-dx/2, xs[-1]+dx/2)

        # Set up plot details
        axes[0, 0].yaxis.tick_right()
        axes[0, 0].set_yscale('log')
        axes[0, 1].set_axis_off()
        axes[0, 2].set_axis_off()
        axes[1, 0].set_ylabel('li')
        axes[1, 1].yaxis.tick_right()
        axes[1, 1].set_yscale('log')
        axes[1, 2].set_axis_off()
        axes[2, 0].set_xlabel('alpha')
        axes[2, 0].set_ylabel('si')
        axes[2, 0].get_shared_y_axes().join(axes[2, 0], axes[2, 1])
        axes[2, 1].set_xlabel('li')
        axes[2, 1].set_yticks([])
        axes[2, 2].set_xlabel('si')
        axes[2, 2].set_yscale('log')
        axes[2, 2].yaxis.tick_right()

        plt.tight_layout()
        plt.subplots_adjust(wspace=0.05, hspace=0.05)

        # # Put colorbar in right position
        cbaxes = inset_axes(axes[0, 2], width="100%", height="3%", loc='center')
        clb = plt.colorbar(im, cax=cbaxes, orientation='horizontal')
        clb.set_label('p-value')

        # Save to subdirectory
        path_to_save = rel_path('./plots/mc/')
        if not os.path.isdir(path_to_save):
            os.mkdir(path_to_save)
        path_to_save += f'{ks_type[3:]}_{survey_name}.pdf'
        plt.savefig(path_to_save)
        plt.clf()

    def make_mesh(self, x_par, y_par, df, ks_type):
        """Make a grid of ks test value for plotting."""
        x_vals = np.sort(df[x_par].unique())
        y_vals = np.sort(df[y_par].unique())

        # To get pixels to line up with labels
        dx = np.diff(x_vals)[0]
        dy = np.diff(y_vals)[0]
        x_vals = np.linspace(min(x_vals), max(x_vals) + dx, len(x_vals)+1)
        y_vals = np.linspace(min(y_vals), max(y_vals) + dy, len(y_vals)+1)
        v = np.zeros([len(x_vals), len(y_vals)])
        pprint(f'{x_par}, {y_par}')

        for i, x_val in enumerate(x_vals):
            for j, y_val in enumerate(y_vals):
                mask = ((df[x_par] == x_val) & (df[y_par] == y_val))
                v[i, j] = df[mask][ks_type].prod(skipna=True, min_count=1)

        return x_vals-dx/2., y_vals-dy/2., v.T

    def make_hist(self, par, df, ks_type):
        """Construct histogram details."""
        par_vals = np.sort(df[par].unique())
        probs = []
        for par_val in par_vals:
            mask = (df[par] == par_val)
            prob = df[mask][ks_type].prod(skipna=True, min_count=1)
            probs.append(prob)

        # Add edges to histograms
        probs = np.array(probs)
        bin_dif = np.diff(par_vals)[-1]
        par_vals = np.insert(par_vals, 0, par_vals[0] - bin_dif)
        par_vals = np.insert(par_vals, len(par_vals), par_vals[-1] + bin_dif)
        probs = np.insert(probs, 0, 0)
        probs = np.insert(probs, len(probs), 0)

        return par_vals, probs

    def apply_test(self, df):
        for c in ('ks_dm', 'ks_snr', 'ks_rate'):
            # Add test data
            df[c] = df[c].apply(lambda v: np.random.random())
        return df

if __name__ == '__main__':
    mc = MonteCarlo(ALPHAS, LIS, SIS, SURVEY_NAMES, POP_SIZE)
    if MAKE:
        mc.gen_runs()
    else:
        mc.get_runs()
    if PLOT:
        mc.plot()
