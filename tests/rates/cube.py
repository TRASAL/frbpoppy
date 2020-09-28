"""Calculate rate cube over alpha, spectral & luminosity index."""
from glob import glob
from matplotlib.widgets import Slider, RadioButtons, CheckButtons
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

from frbpoppy import paths, unpickle, hist, poisson_interval

from alpha_real import EXPECTED

CSV_PATH = os.path.join(paths.populations(), 'cube_rates.csv')

# Set parameters needed for generating a rate cube
GENERATE = False
PLOT = True
SIZE = 1e3
ALPHAS = np.linspace(-1, -2, 5)[1:4]
LIS = np.linspace(-2, 0, 5)  # Luminosity function index
SIS = np.linspace(-2, 2, 5)  # Spectral index
SURVEY_NAMES = ('askap-fly', 'fast-crafts', 'parkes-htru', 'wrst-apertif',
                'arecibo-palfa')

# # Current detection rates
# EXPECTED = {'parkes-htru': [9, 1549 / 0.551 / 24],  # N_frbs, N_days
#             'wrst-apertif': [9, 1100/24],  # 1100 hours
#             'askap-fly': [20, 32840 / 8 / 24],
#             'arecibo-palfa': [1, 24.1],
#             'guppi': [0.4, 81],  # 0.4 is my own assumption
#             'fast-crafts': [1, 1500/24]
#             }


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
                    surv_pop.save()

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

        temp_path = ('./temp.mmap')

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

    df.to_csv(CSV_PATH)


def get_pops(alpha='*', li='*', si='*', survey='*', get_range=False):
    filename = f'complex_alpha_{alpha}_lum_{li}_si_{si}_{survey}.p'
    filter = os.path.join(paths.populations(), filename)
    pop_paths = glob(filter)

    pops = []
    for path in pop_paths:
        if '_for_plotting' not in path:
            pops.append(unpickle(path))
    return pops


def plot():
    # Unpack rates
    df = pd.read_csv(CSV_PATH, index_col=0)
    parameters = ('alpha', 'li', 'si')
    surveys = [c for c in list(df.columns.values) if c not in parameters]
    n_frbs = {k: EXPECTED[k][0] for k in EXPECTED}
    n_days = {k: EXPECTED[k][1] for k in EXPECTED}
    rates = {k: EXPECTED[k][0]/EXPECTED[k][1] for k in EXPECTED}

    # Start plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15*1.5, 9))
    plt.subplots_adjust(right=9/15)

    # Fluence plot
    ax1.set_xlabel('SNR')
    ax1.set_xscale('log')
    ax1.set_ylabel('N > SNR')
    ax1.set_yscale('log')
    ax1.set_xlim(1e0, 1e5)
    ax1.set_ylim(1e-5, 1e3)

    # Rate plot
    ax2.margins(x=0)
    ax2.set_yscale('log')

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # Set up initial plot lines
    # -------------------------
    flnc_lines = []
    rate_lines = []
    expects = []
    intervals = []
    for i, survey in enumerate(surveys):
        # Fluences
        f_line, = ax1.step([0], [0], label=survey, color=colors[i], where='mid')
        flnc_lines.append(f_line)

        # Rates
        r_line, = ax2.plot([0], [0], label=survey, color=colors[i], marker='x')
        rate_lines.append(r_line)

        # # Set up expectation intervals
        interval, = ax2.fill([0, 0, 0, 0], [0, 0, 0, 0], alpha=0.25)
        intervals.append(interval)
        expect, = ax2.plot([0, 0], [0, 0], color=colors[i], linestyle='dashed')
        expects.append(expect)

    def update(*args):
        """Update plot if widgets change."""
        # Which parameter is on the x axis of the rate plot?
        x_axis = x_para.value_selected
        ax2.set_xlabel(x_axis)

        # What values do the parameters have?
        vars = {s.label.get_text(): s.val for s in sliders}

        # Update fluence plot
        for i, line in enumerate(flnc_lines):
            survey = line.get_label()
            pop = get_pops(**vars, survey=survey)[0]

            snr = pop.frbs.snr
            try:
                bins, values = hist(snr, bin_type='log', norm=None)
            except ValueError:
                bins, values = np.array([np.nan]), np.array([np.nan])
            # Cumulative sum
            values = np.cumsum(values[::-1])[::-1]
            # Normalise to area on sky
            values = values * pop.source_rate.f_area
            line.set_xdata(bins)
            line.set_ydata(values)

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

        # Update rate values
        for i, line in enumerate(rate_lines):
            survey = line.get_label()
            x = df[x_axis].unique()
            y = df[mask][survey] / y_scaling

            line.set_xdata(x)
            line.set_ydata(y)

            # Set up expectation intervals
            middle = rates[survey] / real_scaling
            top, bot = poisson_interval(n_frbs[survey], sigma=2)
            top = (top / n_days[survey]) / real_scaling
            bot = (bot / n_days[survey]) / real_scaling
            left = min(x)
            right = max(x)
            xy = list(zip([left, right, right, left], [top, top, bot, bot]))
            intervals[i].set_xy(xy)
            expects[i].set_xdata([left, right])
            expects[i].set_ydata([middle]*2)

        # Fluence plot

        # Rates plot
        ax2.set_xlim(min(x), max(x))
        ax2.set_ylim(1e-2, 1e2)
        ax2.autoscale_view()
        fig.canvas.draw_idle()

    # Ability to select the x-axis
    # ----------------------------
    x_para_pos = plt.axes([9/15+0.03, 1-0.25, 0.15, 0.15])
    x_para = RadioButtons(x_para_pos, parameters, active=0)
    ax2.set_xlabel(parameters[0])
    x_para.on_clicked(update)

    # Make sliders to select values
    # -----------------------------
    offset = 0.4
    sliders = []
    for par in parameters:
        pos = plt.axes([9/15+0.04, offset, 0.15, 0.05])
        step = np.diff(df[par].unique())[0]
        valinit = np.take(df[par], df[par].size//2)
        slider = Slider(pos, par, min(df[par]), max(df[par]),
                        valinit=valinit, valstep=step)
        sliders.append(slider)
        offset -= 0.1

    for slider in sliders:
        slider.on_changed(update)

    # Choose which survey to scale
    # ----------------------------
    scaling_pos = plt.axes([9/15+0.03+0.15, 1-0.25, 0.15, 0.15])
    scaling = RadioButtons(scaling_pos, surveys, active=len(surveys)//2)
    scaling.on_clicked(update)

    # Adapt visibility of lines
    # ----------------------------------------------------------------
    surv_vis_pos = plt.axes([9/15+0.03, 1-0.25*2, 0.15, 0.15])
    visibility = [line.get_visible() for line in rate_lines]
    surv_viz = CheckButtons(surv_vis_pos, surveys, visibility)

    def set_vis(survey):
        index = surveys.index(survey)
        rate_lines[index].set_visible(not rate_lines[index].get_visible())
        flnc_lines[index].set_visible(not flnc_lines[index].get_visible())
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
    ax2.legend(loc='lower right')
    plt.show()


if __name__ == '__main__':
    if GENERATE:
        generate(parallel=True)
    if PLOT:
        plot()
