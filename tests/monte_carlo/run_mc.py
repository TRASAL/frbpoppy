"""Run a Monte Carlo determining best fit parameters.

TODO: Work in progress."""
from frbpoppy import CosmicPopulation, Survey, SurveyPopulation
from frbpoppy import pprint, unpickle, TNS, poisson_interval
from joblib import Parallel, delayed
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import ks_2samp
from scipy.optimize import curve_fit
from tqdm import tqdm
import frbpoppy.paths
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import warnings

from tests.rates.alpha_real import EXPECTED
from tests.convenience import plot_aa_style, rel_path

SURVEY_NAMES = ['parkes-htru', 'chime-frb', 'askap-incoh', 'wsrt-apertif']
NORM_SURV = 'parkes-htru'
ALPHAS = np.linspace(-2, -0.5, 11)
LIS = np.linspace(-2, 0, 11)
SIS = np.linspace(-2, 2, 11)
POP_SIZE = 1e6
MAKE = False
TEST_DATA = False
PLOT = True


class Run:
    """Hold information relevant to a single Monte Carlo run."""

    def __init__(self, alpha, li, si, survey_name=None):
        self.li = li
        self.si = si
        self.alpha = alpha
        self.survey_name = survey_name
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
        self.norm_surv = NORM_SURV
        self.set_up_surveys()
        self.set_up_dir()
        self.get_tns()
        self.runs = []

    def get_tns(self):
        # Only get one-offs
        self.tns = TNS(repeaters=False, mute=True).df

    def set_up_surveys(self):
        """Set up surveys."""
        self.surveys = []
        for name in self.survey_names:
            survey = Survey(name=name)
            survey.set_beam(model='airy', n_sidelobes=1)
            if name in ('chime-frb', 'wsrt-apertif', 'parkes-htru'):
                survey.set_beam(model=name)
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
            n_frbs, n_days = EXPECTED[run.survey_name]
        else:
            n_frbs, n_days = [np.nan, np.nan]

        # Determine ratio of detection rates
        surv_sim_rate = run.rate
        surv_real_rate = n_frbs/n_days

        # Get normalisation properties
        norm_real_n_frbs, norm_real_n_days = EXPECTED[self.norm_surv]
        f = f'mc/complex_alpha_{alpha}_lum_{li}_si_{si}_{self.norm_surv}'
        norm_pop = unpickle(f)
        norm_sim_n_frbs = norm_pop.source_rate.det
        norm_sim_n_days = norm_pop.source_rate.days
        norm_sim_rate = norm_sim_n_frbs / norm_sim_n_days
        norm_real_rate = norm_real_n_frbs / norm_real_n_days

        sim_ratio = surv_sim_rate / norm_sim_rate
        real_ratio = surv_real_rate / norm_real_rate

        # TODO: Update to include information on the Poisson intervals
        run.ks_rate = 1/np.abs(sim_ratio - real_ratio)
        if run.survey_name == self.norm_surv:
            run.ks_rate = np.nan

        mask = (self.tns.survey == run.survey_name)
        run.ks_dm = ks_2samp(surv_pop.frbs.dm, self.tns[mask].dm)[1]
        run.ks_snr = ks_2samp(surv_pop.frbs.snr, self.tns[mask].snr)[1]

        return run

    def plot(self):
        # Convert to a pandas dataframe
        pprint('Plotting')

        df = pd.DataFrame([r.to_dict() for r in self.runs])

        # Suppress warnings on empty slices
        warnings.simplefilter("ignore", category=RuntimeWarning)

        for ks_type in ('ks_dm', 'ks_snr', 'ks_rate'):
            pprint(f'ks-type: {ks_type}')
            for i, group in df.groupby('survey_name'):
                self.plot_run(ks_type, group)
            # A combined plot per ks type
            self.plot_run(ks_type, df)
        # nd one total combined plot
        ks_cols = ['ks_dm', 'ks_snr']
        df['ks_all'] = df[ks_cols].median(axis='columns', skipna=True)
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
        if survey_name == self.norm_surv and ks_type == 'ks_rate':
            return
        pprint(f' - survey: {survey_name}')

        # Add color grids
        args = {'cmap': 'viridis', 'norm': LogNorm(),
                'vmin': 1e-2, 'vmax': 1e0}
        if ks_type == 'ks_rate':
            # Inversting colorbar as smaller is better!
            args = {'cmap': 'viridis', 'norm': LogNorm(),
                    'vmin': 1e-2, 'vmax': 1e0}
        axes[2, 0].pcolormesh(*self.make_mesh('alpha', 'si', df, ks_type),
                              **args)
        axes[1, 0].pcolormesh(*self.make_mesh('alpha', 'li', df, ks_type),
                              **args)
        im = axes[2, 1].pcolormesh(*self.make_mesh('li', 'si', df, ks_type),
                                   **args)

        # Add histograms
        self.plot_hist(axes[0, 0], 'alpha', df, ks_type)
        self.plot_hist(axes[1, 1], 'li', df, ks_type)
        self.plot_hist(axes[2, 2], 'si', df, ks_type)

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
        axes[2, 0].set_xlabel(r'$\alpha$')
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
        cbaxes = inset_axes(axes[0, 2], width="100%", height="3%",
                            loc='center')
        clb = plt.colorbar(im, cax=cbaxes, orientation='horizontal')
        if ks_type == 'ks_rate':
            clb.set_label(f'1/diff rate ratio wrt {self.norm_surv}')
        else:
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
        pprint(f'    - {x_par}, {y_par}')

        for i, x_val in enumerate(x_vals):
            for j, y_val in enumerate(y_vals):
                mask = ((df[x_par] == x_val) & (df[y_par] == y_val))
                v[i, j] = df[mask][ks_type].median(skipna=True)

        return x_vals-dx/2., y_vals-dy/2., v.T

    def plot_hist(self, ax, parameter, df, ks_type):
        bin_centres, values = self.make_hist(parameter, df, ks_type)
        ax.step(bin_centres, values, where='mid', zorder=5)
        # Limited as the bin_centres have extra values to look nice
        bins, hist_fit, coeff, fit_err = self.make_fit(bin_centres[1:-1],
                                                       values[1:-1])
        ax.step(bins, hist_fit, where='mid', color='grey', linewidth=1)
        mu, sigma, norm = coeff
        if not np.isnan(mu):
            ax.axvline(x=mu, linewidth=1, color='grey')
            ax.axvline(x=mu-sigma, linewidth=1, color='grey',
                       linestyle='dashed')
            ax.axvline(x=mu+sigma, linewidth=1, color='grey',
                       linestyle='dashed')
            p = parameter
            if parameter == 'alpha':
                p = r'$\alpha$'
            title = fr'{p}=${mu:.1f}\pm{np.abs(sigma):.1f}$'
            title += '\n'
            title += fr'($\sigma={fit_err:.2e}$)'
            ax.set_title(title,
                         fontsize=10)

    def make_hist(self, par, df, ks_type):
        """Construct histogram details."""
        par_vals = np.sort(df[par].unique())
        probs = []
        for par_val in par_vals:
            mask = (df[par] == par_val)
            prob = df[mask][ks_type].median(skipna=True)
            probs.append(prob)

        # Add edges to histograms
        probs = np.array(probs)
        bin_dif = np.diff(par_vals)[-1]
        par_vals = np.insert(par_vals, 0, par_vals[0] - bin_dif)
        par_vals = np.insert(par_vals, len(par_vals), par_vals[-1] + bin_dif)
        probs = np.insert(probs, 0, 0)
        probs = np.insert(probs, len(probs), 0)

        return par_vals, probs

    def make_fit(self, bin_centres, hist):
        """Fit a function to the constructed median histograms."""
        mask = ~(np.isnan(bin_centres) | np.isnan(hist))
        bin_centres = bin_centres[mask]
        hist = hist[mask]

        def gauss(x, *p):
            mu, sigma, norm = p
            return norm*np.exp(-(x-mu)**2/(2.*sigma**2))

        p0 = [0., 1., 1.]
        try:
            coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
            fit_err = np.sum(np.sqrt(np.diag(var_matrix))**2)
        except (TypeError, RuntimeError) as e:
            return np.nan, np.nan, [np.nan, np.nan, np.nan], np.nan

        # Get the fitted curve
        bins = np.linspace(bin_centres[0], bin_centres[-1], 1000)
        hist_fit = gauss(bins, *coeff)

        return bins, hist_fit, coeff, fit_err

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
