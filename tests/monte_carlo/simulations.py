"""Run Monte Carlo simulations."""
from joblib import Parallel, delayed
from frbpoppy import Survey, CosmicPopulation, SurveyPopulation, pprint
from datetime import datetime
from copy import deepcopy

from glob import glob
import frbpoppy.paths
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import uuid

KEEP_CSV = True
SIZE = 1e7
RUNS = [3]
ALPHA = -1.7
SI = 2.0
LI = -1.0
LI_2 = -0.9
L_MIN = 5e41
L_MAX = 1e45
W_MEAN = 1
W_STD = 1
DM_IGM_SLOPE = 950
DM_HOST = 100


class SimulationOverview:
    """Given values, return uid

    Load from file, or make."""

    def __init__(self, load_csv=True):
        p = frbpoppy.paths.populations()
        self.filename = f'{p}mc/simluation_overview.csv'

        if load_csv and os.path.isfile(self.filename):
            self.load()
        else:
            self.df = pd.DataFrame()

    def load(self):
        self.df = pd.read_csv(self.filename, index_col=0)
        self.df = self.df.loc[:, ~self.df.columns.str.contains('^Unnamed')]

    def save(self):
        self.df.to_csv(self.filename)

    def append(self, df):
        self.df = self.df.append(df, ignore_index=True)

    def map_surveys(self, ix, names):
        mapping = dict(zip(ix, names))
        self.df.replace({"survey": mapping}, inplace=True)


class MonteCarlo:

    def __init__(self, runs=[1], load_csv=True):
        self.survey_names = ['parkes-htru',
                             'chime-frb',
                             'askap-incoh',
                             'wsrt-apertif']
        self.pop_size = SIZE
        self.load_csv = load_csv
        self.runs = runs

        self.survey_ix = [i for i in range(len(self.survey_names))]
        self.surveys = self.set_up_surveys()
        self.so = SimulationOverview(load_csv=self.load_csv)
        self.set_up_dirs()

        for n in self.runs:
            self.so.df = self.so.df[(self.so.df.run != n)]
            if n == 1:
                self.run_1()
            if n == 2:
                self.run_2()
            if n == 3:
                self.run_3()
            if n == 4:
                self.run_4()
            if n == 5:
                self.run_1(5)

    def set_up_surveys(self):
        """Set up surveys."""
        surveys = []
        for name in self.survey_names:
            survey = Survey(name=name)
            survey.set_beam(model='airy', n_sidelobes=1)
            if name in ('chime-frb', 'wsrt-apertif', 'parkes-htru'):
                survey.set_beam(model=name)
            surveys.append(survey)
        return surveys

    def set_up_dirs(self, run_number=None):
        """Create subdirectory for saving populations.

        Returns True if directory had to be set up."""
        f = f'{frbpoppy.paths.populations()}mc/'
        if not os.path.isdir(f):
            os.mkdir(f)
            return True

        if run_number:
            f = f'{frbpoppy.paths.populations()}mc/run_{run_number}/'
            if not os.path.isdir(f):
                os.mkdir(f)
                return True

        return False

    def run_1(self, run_number=1, parallel=True):
        alphas = np.linspace(-2, -0.5, 11)
        sis = np.linspace(-2, 2, 11)
        lis = np.linspace(-2, 0, 11)

        # Put all options into a dataframe
        self.so.df = self.so.df[self.so.df.run != run_number]
        opt = np.meshgrid(alphas, sis, lis, self.survey_ix)
        options = np.array(opt).T.reshape(-1, 4)
        df = pd.DataFrame(options, columns=('alpha', 'si', 'li', 'survey'))
        df['run'] = run_number
        df['uuid'] = [uuid.uuid4() for _ in range(len(df.index))]
        df['date'] = datetime.today()
        self.so.append(df)
        self.so.map_surveys(self.survey_ix, self.survey_names)
        self.so.save()

        # Remove previous run of the same number
        if not self.set_up_dirs(run_number=run_number):
            fs = f'{frbpoppy.paths.populations()}mc/run_{run_number}/*'
            for f in glob(fs):
                os.remove(f)

        def iter_alpha(i, parallel=None):
            alpha = alphas[i]
            pop = CosmicPopulation.complex(self.pop_size)

            # # Basically a sanity check
            # if run_number == 5:
            #     pop.set_w(model='lognormal', mean=W_MEAN, std=W_STD)
            #     pop.set_dm_igm(model='ioka', slope=DM_IGM_SLOPE)
            #     pop.set_dm_host(model='constant', value=DM_HOST)

            pop.set_dist(model='vol_co', z_max=1.0, alpha=alpha)
            pop.set_lum(model='constant', value=1)
            pop.generate()

            for si in sis:
                pop.set_si(model='constant', value=si)
                pop.gen_si()

                for li in lis:
                    pop.set_lum(model='powerlaw',
                                low=1e40,
                                high=1e45, power=li)

                    if run_number == 5:
                        pop.set_lum(model='powerlaw', low=L_MIN,
                                    high=L_MAX, index=li)

                    pop.gen_lum()

                    for survey in self.surveys:
                        surv_pop = SurveyPopulation(pop, survey)

                        # Get unique identifier
                        mask = (self.so.df.run == run_number)
                        mask &= (self.so.df.alpha == alpha)
                        mask &= (self.so.df.si == si)
                        mask &= (self.so.df.li == li)
                        mask &= (self.so.df.survey == survey.name)
                        uuid = self.so.df[mask].uuid.iloc[0]
                        surv_pop.name = f'mc/run_{run_number}/{uuid}'
                        surv_pop.save()

        if parallel:
            n_cpu = min([3, os.cpu_count() - 1])
            pprint(f'{os.cpu_count()} CPUs available')
            r = range(len(alphas))
            Parallel(n_jobs=n_cpu)(delayed(iter_alpha)(i) for i in tqdm(r))
        else:
            [iter_alpha(i) for i in tqdm(range(len(alphas)))]

    def run_2(self, run_number=2, parallel=True):
        lis = np.linspace(-1.5, 0, 11)
        lum_mins = 10**np.linspace(34, 45, 11)
        lum_maxs = 10**np.linspace(34, 45, 11)

        # Put all options into a dataframe
        self.so.df = self.so.df[self.so.df.run != 2]
        opt = np.meshgrid(lis, lum_mins, lum_maxs, self.survey_ix)
        options = np.array(opt).T.reshape(-1, 4)
        cols = ('li', 'lum_min', 'lum_max', 'survey')
        df = pd.DataFrame(options, columns=cols)
        df['run'] = run_number
        df['uuid'] = [uuid.uuid4() for _ in range(len(df.index))]
        df['date'] = datetime.today()
        df = df[~(df.lum_max < df.lum_min)]
        self.so.append(df)
        self.so.map_surveys(self.survey_ix, self.survey_names)
        self.so.save()

        # Remove previous run of the same number
        if not self.set_up_dirs(run_number=run_number):
            fs = f'{frbpoppy.paths.populations()}mc/run_{run_number}/*'
            for f in glob(fs):
                os.remove(f)

        pop = CosmicPopulation.complex(self.pop_size)
        pop.set_dist(model='vol_co', z_max=1.0, alpha=ALPHA)
        pop.set_si(model='constant', value=SI)
        pop.set_lum(model='constant', value=1)
        pop.generate()

        def adapt_pop(e):
            li, lum_min, lum_max = e
            if lum_max < lum_min:
                return
            t_pop = deepcopy(pop)
            t_pop.set_lum(model='powerlaw', low=lum_min, high=lum_max,
                          power=li)
            t_pop.gen_lum()

            for survey in self.surveys:
                surv_pop = SurveyPopulation(t_pop, survey)

                # Get unique identifier
                mask = (self.so.df.run == run_number)
                mask &= (self.so.df.li == li)
                mask &= (self.so.df.lum_min == lum_min)
                mask &= (self.so.df.lum_max == lum_max)
                mask &= (self.so.df.survey == survey.name)
                uuid = self.so.df[mask].uuid.iloc[0]
                surv_pop.name = f'mc/run_{run_number}/{uuid}'
                surv_pop.save()

        n_cpu = min([3, os.cpu_count() - 1])
        pprint(f'{os.cpu_count()} CPUs available')
        mg = np.meshgrid(lis, lum_mins, lum_maxs)
        loop_over = np.array(mg).T.reshape(-1, 3)
        if parallel:
            Parallel(n_jobs=n_cpu)(delayed(adapt_pop)(e) for e in tqdm(loop_over))
        else:
            [adapt_pop(e) for e in tqdm(loop_over)]

    def run_3(self, run_number=3, parallel=True):
        w_means = 10**np.linspace(-2, 2, 11)  # TODO!
        w_stds = np.linspace(0, 3, 11)  # TODO!

        # Put all options into a dataframe
        self.so.df = self.so.df[self.so.df.run != 3]
        opt = np.meshgrid(w_means, w_stds, self.survey_ix)
        options = np.array(opt).T.reshape(-1, 3)
        cols = ('w_mean', 'w_std', 'survey')
        df = pd.DataFrame(options, columns=cols)
        df['run'] = run_number
        df['uuid'] = [uuid.uuid4() for _ in range(len(df.index))]
        df['date'] = datetime.today()
        self.so.append(df)
        self.so.map_surveys(self.survey_ix, self.survey_names)
        self.so.save()

        # Remove previous run of the same number
        if not self.set_up_dirs(run_number=run_number):
            fs = f'{frbpoppy.paths.populations()}mc/run_{run_number}/*'
            for f in glob(fs):
                os.remove(f)

        pop = CosmicPopulation.complex(self.pop_size)
        pop.set_dist(model='vol_co', z_max=1.0, alpha=ALPHA)
        pop.set_si(model='constant', value=SI)
        pop.set_lum(model='powerlaw', low=L_MIN, high=L_MAX, index=LI_2)
        pop.generate()

        def adapt_pop(e):
            w_mean, w_std = e
            t_pop = deepcopy(pop)
            t_pop.set_w(model='lognormal', mean=w_mean, std=w_std)
            t_pop.gen_w()

            for survey in self.surveys:
                surv_pop = SurveyPopulation(t_pop, survey)

                # Get unique identifier
                mask = (self.so.df.run == run_number)
                mask &= (self.so.df.w_mean == w_mean)
                mask &= (self.so.df.w_std == w_std)
                mask &= (self.so.df.survey == survey.name)
                uuid = self.so.df[mask].uuid.iloc[0]
                surv_pop.name = f'mc/run_{run_number}/{uuid}'
                surv_pop.save()

        n_cpu = min([3, os.cpu_count() - 1])
        pprint(f'{os.cpu_count()} CPUs available')
        mg = np.meshgrid(w_means, w_stds)
        loop_over = np.array(mg).T.reshape(-1, 2)
        if parallel:
            Parallel(n_jobs=n_cpu)(delayed(adapt_pop)(e) for e in tqdm(loop_over))
        else:
            [adapt_pop(e) for e in tqdm(loop_over)]

    def run_4(self, run_number=4, parallel=True):
        dm_igm_slopes = np.linspace(800, 1200, 11)
        dm_hosts = np.linspace(0, 500, 11)

        # Put all options into a dataframe
        self.so.df = self.so.df[self.so.df.run != 4]
        opt = np.meshgrid(dm_igm_slopes, dm_hosts, self.survey_ix)
        options = np.array(opt).T.reshape(-1, 3)
        cols = ('dm_igm_slope', 'dm_host', 'survey')
        df = pd.DataFrame(options, columns=cols)
        df['run'] = run_number
        df['uuid'] = [uuid.uuid4() for _ in range(len(df.index))]
        df['date'] = datetime.today()
        self.so.append(df)
        self.so.map_surveys(self.survey_ix, self.survey_names)
        self.so.save()

        # Remove previous run of the same number
        if not self.set_up_dirs(run_number=run_number):
            fs = f'{frbpoppy.paths.populations()}mc/run_{run_number}/*'
            for f in glob(fs):
                os.remove(f)

        pop = CosmicPopulation.complex(self.pop_size)
        pop.set_dist(model='vol_co', z_max=1.0, alpha=ALPHA)
        pop.set_si(model='constant', value=SI)
        pop.set_lum(model='powerlaw', low=L_MIN, high=L_MAX, index=LI_2)
        pop.set_w(model='lognormal', mean=W_MEAN, std=W_STD)
        pop.generate()

        def adapt_pop(e):
            dm_igm_slope, dm_host = e
            t_pop = deepcopy(pop)
            t_pop.set_dm_igm(model='ioka', slope=dm_igm_slope)
            t_pop.gen_dm_igm()
            t_pop.set_dm_host(model='constant', value=dm_host)
            t_pop.gen_dm_host()
            t_pop.frbs.dm = t_pop.frbs.dm_mw + t_pop.frbs.dm_igm + t_pop.frbs.dm_host

            for survey in self.surveys:
                surv_pop = SurveyPopulation(t_pop, survey)

                # Get unique identifier
                mask = (self.so.df.run == run_number)
                mask &= (self.so.df.dm_igm_slope == dm_igm_slope)
                mask &= (self.so.df.dm_host == dm_host)
                mask &= (self.so.df.survey == survey.name)
                uuid = self.so.df[mask].uuid.iloc[0]
                surv_pop.name = f'mc/run_{run_number}/{uuid}'
                surv_pop.save()

        n_cpu = min([3, os.cpu_count() - 1])
        pprint(f'{os.cpu_count()} CPUs available')
        mg = np.meshgrid(dm_igm_slopes, dm_hosts)
        loop_over = np.array(mg).T.reshape(-1, 2)
        if parallel:
            Parallel(n_jobs=n_cpu)(delayed(adapt_pop)(e) for e in tqdm(loop_over))
        else:
            [adapt_pop(e) for e in tqdm(loop_over)]


if __name__ == '__main__':
    MonteCarlo(load_csv=KEEP_CSV,
               runs=RUNS)
