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

POP_SIZE = 5e7


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

    def __init__(self, pop_size=1e2, load_csv=True):
        self.survey_names = ['parkes-htru',
                             'chime-frb',
                             'askap-incoh',
                             'wsrt-apertif']
        self.load_csv = load_csv
        self.pop_size = pop_size
        self.survey_ix = [i for i in range(len(self.survey_names))]
        self.surveys = self.set_up_surveys()
        self.so = SimulationOverview(load_csv=self.load_csv)
        self.set_up_dirs()

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

    def set_up_dirs(self, run=np.nan):
        """Create subdirectory for saving populations.

        Returns True if directory had to be set up."""
        f = f'{frbpoppy.paths.populations()}mc/'
        if not os.path.isdir(f):
            os.mkdir(f)
            return True

        if not np.isnan(run):
            f = f'{frbpoppy.paths.populations()}mc/run_{run}/'
            if not os.path.isdir(f):
                os.mkdir(f)
                return True

        return False

    def gen_par_set_1(self,
                      parallel=True,
                      lum_min=np.nan,
                      lum_max=np.nan,
                      w_mean=np.nan,
                      w_std=np.nan,
                      dm_igm_slope=np.nan,
                      dm_host=np.nan,
                      run=0):
        alphas = np.linspace(-2.5, -1, 11)
        sis = np.linspace(-2, 2, 11)
        lis = np.linspace(-2, 0, 11)

        # Put all options into a dataframe
        if 'run' in self.so.df:
            self.so.df = self.so.df[self.so.df.run != run]
        opt = np.meshgrid(alphas, sis, lis, self.survey_ix)
        options = np.array(opt).T.reshape(-1, 4)
        df = pd.DataFrame(options, columns=('alpha', 'si', 'li', 'survey'))
        df['run'] = run
        df['par_set'] = 1
        df['uuid'] = [uuid.uuid4() for _ in range(len(df.index))]
        df['date'] = datetime.today()
        self.so.append(df)
        self.so.map_surveys(self.survey_ix, self.survey_names)
        self.so.save()

        # Remove previous par_set of the same number
        if not self.set_up_dirs(run=run):
            fs = f'{frbpoppy.paths.populations()}mc/run_{run}/*'
            for f in glob(fs):
                os.remove(f)

        def iter_alpha(i):
            alpha = alphas[i]
            pop = CosmicPopulation.complex(self.pop_size)

            pop.set_dist(model='vol_co', z_max=1.0, alpha=alpha)
            pop.set_lum(model='constant', value=1)

            if not np.isnan(w_mean):
                pop.set_w(model='lognormal', mean=w_mean, std=w_std)
            if not np.isnan(dm_igm_slope):
                pop.set_dm_igm(model='ioka', slope=dm_igm_slope)
                pop.set_dm_host(model='constant', value=dm_host)

            pop.generate()

            for si in sis:
                pop.set_si(model='constant', value=si)
                pop.gen_si()

                for li in lis:
                    pop.set_lum(model='powerlaw',
                                low=1e40,
                                high=1e45, power=li)

                    if not np.isnan(lum_min):
                        pop.set_lum(model='powerlaw', low=lum_min,
                                    high=lum_max, index=li)

                    pop.gen_lum()

                    for survey in self.surveys:
                        surv_pop = SurveyPopulation(pop, survey)

                        # Get unique identifier
                        mask = (self.so.df.par_set == 1)
                        mask &= (self.so.df.run == run)
                        mask &= (self.so.df.alpha == alpha)
                        mask &= (self.so.df.si == si)
                        mask &= (self.so.df.li == li)
                        mask &= (self.so.df.survey == survey.name)
                        uuid = self.so.df[mask].uuid.iloc[0]
                        surv_pop.name = f'mc/run_{run}/{uuid}'
                        surv_pop.save()

        if parallel:
            n_cpu = min([3, os.cpu_count() - 1])
            pprint(f'{os.cpu_count()} CPUs available')
            r = range(len(alphas))
            Parallel(n_jobs=n_cpu)(delayed(iter_alpha)(i) for i in tqdm(r))
        else:
            [iter_alpha(i) for i in tqdm(range(len(alphas)))]

    def gen_par_set_2(self,
                      parallel=True,
                      alpha=-1.5,
                      si=0,
                      w_mean=np.nan,
                      w_std=np.nan,
                      dm_igm_slope=np.nan,
                      dm_host=np.nan,
                      run=np.nan):
        lis = np.linspace(-1.5, 0, 11)
        lum_mins = 10**np.linspace(38, 46, 11)
        lum_maxs = 10**np.linspace(38, 46, 11)

        # Put all options into a dataframe
        self.so.df = self.so.df[self.so.df.run != run]
        opt = np.meshgrid(lis, lum_mins, lum_maxs, self.survey_ix)
        options = np.array(opt).T.reshape(-1, 4)
        cols = ('li', 'lum_min', 'lum_max', 'survey')
        df = pd.DataFrame(options, columns=cols)
        df['par_set'] = 2
        df['run'] = run
        df['uuid'] = [uuid.uuid4() for _ in range(len(df.index))]
        df['date'] = datetime.today()
        df = df[~(df.lum_max < df.lum_min)]
        self.so.append(df)
        self.so.map_surveys(self.survey_ix, self.survey_names)
        self.so.save()

        # Remove previous par_set of the same number
        if not self.set_up_dirs(run=run):
            fs = f'{frbpoppy.paths.populations()}mc/run_{run}/*'
            for f in glob(fs):
                os.remove(f)

        pop = CosmicPopulation.complex(self.pop_size)
        if not np.isnan(alpha):
            pop.set_dist(model='vol_co', z_max=1.0, alpha=alpha)
            pop.set_si(model='constant', value=si)
            pop.set_lum(model='constant', value=1)
        if not np.isnan(w_mean):
            pop.set_w(model='lognormal', mean=w_mean, std=w_std)
        if not np.isnan(dm_igm_slope):
            pop.set_dm_igm(model='ioka', slope=dm_igm_slope)
            pop.set_dm_host(model='constant', value=dm_host)

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
                mask = (self.so.df.par_set == 2)
                mask &= (self.so.df.run == run)
                mask &= (self.so.df.li == li)
                mask &= (self.so.df.lum_min == lum_min)
                mask &= (self.so.df.lum_max == lum_max)
                mask &= (self.so.df.survey == survey.name)
                uuid = self.so.df[mask].uuid.iloc[0]
                surv_pop.name = f'mc/run_{run}/{uuid}'
                surv_pop.save()

        n_cpu = min([3, os.cpu_count() - 1])
        pprint(f'{os.cpu_count()} CPUs available')
        mg = np.meshgrid(lis, lum_mins, lum_maxs)
        loop = np.array(mg).T.reshape(-1, 3)
        if parallel:
            Parallel(n_jobs=n_cpu)(delayed(adapt_pop)(e) for e in tqdm(loop))
        else:
            [adapt_pop(e) for e in tqdm(loop)]

    def gen_par_set_3(self,
                      parallel=True,
                      alpha=-1.5,
                      si=0,
                      li=-1,
                      lum_min=1e40,
                      lum_max=1e40,
                      dm_igm_slope=np.nan,
                      dm_host=np.nan,
                      run=np.nan):
        w_means = 10**np.linspace(-3, 1, 11)
        w_stds = np.linspace(0, 3, 11)

        # Put all options into a dataframe
        self.so.df = self.so.df[self.so.df.run != run]
        opt = np.meshgrid(w_means, w_stds, self.survey_ix)
        options = np.array(opt).T.reshape(-1, 3)
        cols = ('w_mean', 'w_std', 'survey')
        df = pd.DataFrame(options, columns=cols)
        df['run'] = run
        df['par_set'] = 3
        df['uuid'] = [uuid.uuid4() for _ in range(len(df.index))]
        df['date'] = datetime.today()
        self.so.append(df)
        self.so.map_surveys(self.survey_ix, self.survey_names)
        self.so.save()

        # Remove previous par_set of the same number
        if not self.set_up_dirs(run=run):
            fs = f'{frbpoppy.paths.populations()}mc/run_{run}/*'
            for f in glob(fs):
                os.remove(f)

        pop = CosmicPopulation.complex(self.pop_size)

        if not np.isnan(alpha):
            pop.set_dist(model='vol_co', z_max=1.0, alpha=alpha)
            pop.set_si(model='constant', value=si)
        if not np.isnan(lum_min):
            pop.set_lum(model='powerlaw', low=lum_min, high=lum_max, index=li)
        if not np.isnan(dm_igm_slope):
            pop.set_dm_igm(model='ioka', slope=dm_igm_slope)
            pop.set_dm_host(model='constant', value=dm_host)

        pop.generate()

        def adapt_pop(e):
            w_mean, w_std = e
            t_pop = deepcopy(pop)
            t_pop.set_w(model='lognormal', mean=w_mean, std=w_std)
            t_pop.gen_w()

            for survey in self.surveys:
                surv_pop = SurveyPopulation(t_pop, survey)

                # Get unique identifier
                mask = (self.so.df.par_set == 3)
                mask &= (self.so.df.run == run)
                mask &= (self.so.df.run == run)
                mask &= (self.so.df.w_mean == w_mean)
                mask &= (self.so.df.w_std == w_std)
                mask &= (self.so.df.survey == survey.name)
                uuid = self.so.df[mask].uuid.iloc[0]
                surv_pop.name = f'mc/run_{run}/{uuid}'
                surv_pop.save()

        n_cpu = min([3, os.cpu_count() - 1])
        pprint(f'{os.cpu_count()} CPUs available')
        mg = np.meshgrid(w_means, w_stds)
        loop = np.array(mg).T.reshape(-1, 2)
        if parallel:
            Parallel(n_jobs=n_cpu)(delayed(adapt_pop)(e) for e in tqdm(loop))
        else:
            [adapt_pop(e) for e in tqdm(loop)]

    def gen_par_set_4(self,
                      parallel=True,
                      alpha=-1.5,
                      si=0,
                      li=-1,
                      lum_min=1e40,
                      lum_max=1e40,
                      w_mean=np.nan,
                      w_std=np.nan,
                      run=np.nan):
        dm_igm_slopes = np.linspace(800, 1200, 11)
        dm_hosts = np.linspace(0, 500, 11)

        # Put all options into a dataframe
        self.so.df = self.so.df[self.so.df.run != run]
        opt = np.meshgrid(dm_igm_slopes, dm_hosts, self.survey_ix)
        options = np.array(opt).T.reshape(-1, 3)
        cols = ('dm_igm_slope', 'dm_host', 'survey')
        df = pd.DataFrame(options, columns=cols)
        df['run'] = run
        df['par_set'] = 4
        df['uuid'] = [uuid.uuid4() for _ in range(len(df.index))]
        df['date'] = datetime.today()
        self.so.append(df)
        self.so.map_surveys(self.survey_ix, self.survey_names)
        self.so.save()

        # Remove previous par_set of the same number
        if not self.set_up_dirs(run=run):
            fs = f'{frbpoppy.paths.populations()}mc/run_{run}/*'
            for f in glob(fs):
                os.remove(f)

        pop = CosmicPopulation.complex(self.pop_size)

        if not np.isnan(alpha):
            pop.set_dist(model='vol_co', z_max=1.0, alpha=alpha)
            pop.set_si(model='constant', value=si)
        if not np.isnan(lum_min):
            pop.set_lum(model='powerlaw', low=lum_min, high=lum_max, index=li)
        if not np.isnan(w_mean):
            pop.set_w(model='lognormal', mean=w_mean, std=w_std)
        pop.generate()

        def adapt_pop(e):
            dm_igm_slope, dm_host = e
            t_pop = deepcopy(pop)
            t_pop.set_dm_igm(model='ioka', slope=dm_igm_slope)
            t_pop.gen_dm_igm()
            t_pop.set_dm_host(model='constant', value=dm_host)
            t_pop.gen_dm_host()
            t_pop.frbs.dm = t_pop.frbs.dm_mw + t_pop.frbs.dm_igm
            t_pop.frbs.dm += t_pop.frbs.dm_host

            for survey in self.surveys:
                surv_pop = SurveyPopulation(t_pop, survey)

                # Get unique identifier
                mask = (self.so.df.par_set == 4)
                mask &= (self.so.df.run == run)
                mask &= (self.so.df.dm_igm_slope == dm_igm_slope)
                mask &= (self.so.df.dm_host == dm_host)
                mask &= (self.so.df.survey == survey.name)
                uuid = self.so.df[mask].uuid.iloc[0]
                surv_pop.name = f'mc/run_{run}/{uuid}'
                surv_pop.save()

        n_cpu = min([4, os.cpu_count() - 1])
        pprint(f'{os.cpu_count()} CPUs available')
        mg = np.meshgrid(dm_igm_slopes, dm_hosts)
        loop = np.array(mg).T.reshape(-1, 2)
        if parallel:
            Parallel(n_jobs=n_cpu)(delayed(adapt_pop)(e) for e in tqdm(loop))
        else:
            [adapt_pop(e) for e in tqdm(loop)]
