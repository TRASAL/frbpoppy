"""Calculate rate cube over alpha, spectral & luminosity index."""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons, CheckButtons
import numpy as np
import os
from scipy.stats import chi2, norm
from scipy.integrate import quad

from tests.convenience import rel_path

GENERATE = False
SIZE = 1e8
ALPHAS = np.array([-1, -1.5, -2])
LIS = np.array([-2, -1, 0])  # Luminosity function index
SIS = np.array([-4, 0, 4])  # Spectral index

SURVEY_NAMES = ('askap-fly', 'fast', 'htru', 'apertif', 'palfa')

EXPECTED = {'htru': [9, 0.551 / 1549 / 24],  # N_frbs, N_days
            'apertif': [9, 1100/24],  # 1100 hours
            'askap-fly': [20, 32840 * 8 / 24],
            'palfa': [1, 24.1],
            'guppi': [0.4, 81],  # 0.4 is my own assumption
            'fast': [1, 1500/24]
            }


def generate(parallel=False):
    from joblib import Parallel, delayed
    from frbpoppy import CosmicPopulation, Survey, pprint
    from frbpoppy import SurveyPopulation
    from tqdm import tqdm

    # Set up rate dataframe
    paras = [ALPHAS, LIS, SIS]
    vals = np.array(np.meshgrid(*paras)).T.reshape(-1, len(paras))
    cols = ['alpha', 'li', 'si']
    df = pd.DataFrame(vals, columns=cols)

    # Set up surveys
    surveys = []
    for name in SURVEY_NAMES:
        survey = Survey(name=name)
        survey.set_beam(model='airy', n_sidelobes=1)
        surveys.append(survey)
        df[name] = np.nan

    def iter_alpha(i, surveys=surveys, parallel=None):
        alpha = ALPHAS[i]
        pop = CosmicPopulation.complex(SIZE)
        pop.set_dist(model='vol_co', z_max=1.0, alpha=alpha)
        pop.set_lum(model='powerlaw', low=1e40, high=1e45, power=-1)
        pop.generate()

        for li in LIS:
            pop.set_lum(model='powerlaw', low=1e40, high=1e45, power=li)
            pop.gen_lum()

            for si in SIS:
                pop.set_si(model='constant', value=si)
                pop.gen_si()

                pop.name = f'complex_alpha_{alpha}_lum_{li}_si_{si}'

                for survey in surveys:
                    surv_pop = SurveyPopulation(pop, survey)
                    print(surv_pop.name)

                    sr = surv_pop.source_rate
                    rate = sr.det / sr.days
                    mask = (df.alpha == alpha) & (df.li == li) & (df.si == si)

                    if parallel is not None:
                        i = df[mask].index
                        j = SURVEY_NAMES.index(survey.name)
                        parallel[i, j] = rate
                    else:
                        df.loc[mask, survey.name] = rate

    if parallel:
        n_cpu = min([3, os.cpu_count() - 1])
        pprint(f'{os.cpu_count()} CPUs available')
        r = range(len(ALPHAS))

        temp_path = rel_path('./plots/temp.mmap')

        # Make a temp memmap to have a sharedable memory object
        temp = np.memmap(temp_path, dtype=np.float64,
                         shape=(len(vals), len(SURVEY_NAMES)),
                         mode='w+')

        Parallel(n_jobs=n_cpu)(delayed(iter_alpha)(i, parallel=temp) for i in tqdm(r))

        for name in SURVEY_NAMES:
            col = SURVEY_NAMES.index(name)
            df[name] = temp[:, col]
    else:
        for i in tqdm(range(len(ALPHAS)), desc='Alphas'):
            iter_alpha(i)

    df.to_csv(rel_path('plots/cube_rates.csv'))


def poisson_interval(k, sigma=1):
    """
    Use chi-squared info to get the poisson interval.

    Give a number of observed events, which range of observed events would have
    been just as likely given a particular interval?

    Based off https://stackoverflow.com/questions/14813530/
    poisson-confidence-interval-with-numpy
    """
    gauss = norm(0, 1).pdf
    a = 1 - quad(gauss, -sigma, sigma, limit=1000)[0]
    low, high = (chi2.ppf(a/2, 2*k) / 2, chi2.ppf(1-a/2, 2*k + 2) / 2)
    if k == 0:
        low = 0.0

    return low, high


def plot():
    # Unpack rates
    df = pd.read_csv(rel_path('plots/cube_rates.csv'), index_col=0)
    parameters = ('alpha', 'li', 'si')
    surveys = [c for c in list(df.columns.values) if c not in parameters]
    n_frbs = {k: EXPECTED[k][0] for k in EXPECTED}
    n_days = {k: EXPECTED[k][1] for k in EXPECTED}
    rates = {k: EXPECTED[k][0]/EXPECTED[k][1] for k in EXPECTED}

    # Start plot
    fig, ax = plt.subplots(figsize=(15, 9))
    plt.subplots_adjust(right=9/15)
    ax.margins(x=0)
    ax.set_yscale('log')
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # Set up initial plot lines
    # -------------------------
    lines = []
    expects = []
    intervals = []
    for i, survey in enumerate(surveys):
        line, = plt.plot([0], [0], label=survey, color=colors[i], marker='x')
        lines.append(line)

        # # Set up expectation intervals
        interval, = plt.fill([0, 0, 0, 0], [0, 0, 0, 0], alpha=0.25)
        intervals.append(interval)
        expect, = plt.plot([0, 0], [0, 0], color=colors[i], linestyle='dashed')
        expects.append(expect)

    def update(*args):
        """Update plot if widgets change."""
        # Which parameter is on the x axis?
        x_axis = x_para.value_selected
        ax.set_xlabel(x_axis)

        # What values do the parameters have?
        vars = {s.label.get_text(): s.val for s in sliders}

        # Figure out new condition
        masks = []
        for p in parameters:
            if p != x_axis:
                masks.append((df[p] == vars[p]))
        mask = np.bitwise_and.reduce(np.array(masks))

        # Scale to which survey
        scale_survey = scaling.value_selected
        y_scaling = df[mask][scale_survey].values
        real_scaling = rates[scale_survey]

        for i, line in enumerate(lines):
            survey = line.get_label()
            x = df[x_axis].unique()
            y = df[mask][survey] / y_scaling

            line.set_xdata(x)
            line.set_ydata(y)

            # Set up expectation intervals
            middle = rates[survey] / real_scaling
            top, bot = poisson_interval(n_frbs[survey])
            top = (top / n_days[survey]) / real_scaling
            bot = (bot / n_days[survey]) / real_scaling
            left = min(x)
            right = max(x)
            xy = list(zip([left, right, right, left], [top, top, bot, bot]))
            intervals[i].set_xy(xy)
            expects[i].set_xdata([left, right])
            expects[i].set_ydata([middle]*2)

        # ax.relim()
        ax.set_xlim(min(x), max(x))
        ax.set_ylim(1e-2, 1e2)
        ax.autoscale_view()
        fig.canvas.draw_idle()

    # Ability to select the x-axis
    # ----------------------------
    x_para_pos = plt.axes([9/15+0.03, 1-0.25, 0.15, 0.15])
    x_para = RadioButtons(x_para_pos, parameters, active=0)
    ax.set_xlabel(parameters[0])
    x_para.on_clicked(update)

    # Make sliders to select values
    # -----------------------------
    offset = 0.4
    sliders = []
    for par in parameters:
        pos = plt.axes([9/15+0.04, offset, 0.15, 0.05])
        step = np.diff(df[par].unique())[0]
        slider = Slider(pos, par, min(df[par]), max(df[par]),
                        valinit=df[par].unique()[1], valstep=step)
        sliders.append(slider)
        offset -= 0.1

    for slider in sliders:
        slider.on_changed(update)

    # Choose which survey to scale
    # ----------------------------
    scaling_pos = plt.axes([9/15+0.03+0.15, 1-0.25, 0.15, 0.15])
    scaling = RadioButtons(scaling_pos, surveys, active=len(surveys)-1)
    scaling.on_clicked(update)

    # Adapt visibility of lines
    # ----------------------------------------------------------------
    surv_vis_pos = plt.axes([9/15+0.03, 1-0.25*2, 0.15, 0.15])
    visibility = [line.get_visible() for line in lines]
    surv_viz = CheckButtons(surv_vis_pos, surveys, visibility)

    def set_vis(survey):
        index = surveys.index(survey)
        lines[index].set_visible(not lines[index].get_visible())
        plt.draw()

    surv_viz.on_clicked(set_vis)

    # Adapt visibility of plotted intervals
    # ----------------------------------------------------------------
    interval_vis_pos = plt.axes([9/15+0.03+0.15, 1-0.25*2, 0.15, 0.15])
    visibility = [interval.get_visible() for interval in intervals]
    interval_viz = CheckButtons(interval_vis_pos, surveys, visibility)

    def set_interval(survey):
        index = surveys.index(survey)
        intervals[index].set_visible(not intervals[index].get_visible())
        expects[index].set_visible(not expects[index].get_visible())
        plt.draw()

    interval_viz.on_clicked(set_interval)

    update()
    ax.legend(loc='lower right')
    plt.show()


if __name__ == '__main__':
    if GENERATE:
        generate(parallel=True)

    plot()
