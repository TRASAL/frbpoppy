"""Run Monte Carlo simulations."""
from joblib import Parallel, delayed
from frbpoppy import Survey, CosmicPopulation, SurveyPopulation, pprint

from glob import glob
import frbpoppy.paths
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import uuid


class SimulationOverview:
    """Given values, return uid

    Load from file, or make."""

    def __init__(self, load_csv=False):
        p = frbpoppy.paths.populations()
        self.filename = f'{p}mc/simluation_overview.csv'

        if load_csv:
            self.load()
        else:
            self.df = pd.DataFrame()

    def load(self):
        self.df = pd.read_csv(self.filename)

    def save(self):
        self.df.to_csv(self.filename)

    def append(self, df):
        self.df = self.df.append(df)

    def map_surveys(self, ix, names):
        mapping = dict(zip(ix, names))
        self.df.replace({"survey": mapping}, inplace=True)


class MonteCarlo:

    def __init__(self, runs=[1], load_csv=True):
        self.survey_names = ['parkes-htru',
                             'chime-frb',
                             'askap-incoh',
                             'wsrt-apertif']
        self.pop_size = 1e7
        self.load_csv = load_csv
        self.runs = runs

        self.survey_ix = [i for i in range(len(self.survey_names))]
        self.surveys = self.set_up_surveys()
        self.so = SimulationOverview(load_csv=self.load_csv)
        self.set_up_dir()

        if 1 in self.runs:
            self.run_1()

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

    def set_up_dir(self):
        """Create subdirectory for saving populations."""
        f = f'{frbpoppy.paths.populations()}mc/'
        if not os.path.isdir(f):
            os.mkdir(f)

    def run_1(self, run_number=1, parallel=False):
        alphas = np.linspace(-2, -0.5, 11)
        sis = np.linspace(-2, 2, 11)
        lis = np.linspace(-2, 0, 11)

        # Put all options into a dataframe
        opt = np.meshgrid(alphas, sis, lis, self.survey_ix)
        options = np.array(opt).T.reshape(-1, 4)
        df = pd.DataFrame(options, columns=('alpha', 'si', 'li', 'survey'))
        df['run'] = run_number
        df['uuid'] = [uuid.uuid4() for _ in range(len(df.index))]
        self.so.append(df)
        self.so.map_surveys(self.survey_ix, self.survey_names)
        self.so.save()

        # Remove previous run of the same number
        fs = f'{frbpoppy.paths.populations()}mc/run_{run_number}_*'
        for f in glob(fs):
            os.remove(f)

        def iter_alpha(i, parallel=None):
            alpha = alphas[i]
            pop = CosmicPopulation.complex(self.pop_size)
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
                        surv_pop.name = f'mc/run_{run_number}_{uuid}'
                        surv_pop.save()

        if parallel:
            n_cpu = min([3, os.cpu_count() - 1])
            pprint(f'{os.cpu_count()} CPUs available')
            r = range(len(alphas))
            Parallel(n_jobs=n_cpu)(delayed(iter_alpha)(i) for i in tqdm(r))
        else:
            [iter_alpha(i) for i in tqdm(range(len(alphas)))]

    def run_2():
        pass

    def run_3():
        pass

    def run_4():
        pass


if __name__ == '__main__':
    MonteCarlo(load_csv=False,
               runs=[1])
