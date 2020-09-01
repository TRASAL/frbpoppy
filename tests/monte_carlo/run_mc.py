"""Run a Monte Carlo determining best fit parameters."""
from frbpoppy import CosmicPopulation, Survey, pprint
from frbpoppy import SurveyPopulation
from joblib import Parallel, delayed
from tqdm import tqdm
import numpy as np
import pandas as pd
import os

SURVEY_NAMES = ['apertif', 'askap-fly', 'htru', 'chime']
ALPHAS = np.linspace(-2, -0.5, 6)
LIS = np.linspace(-2, 0, 5)
SIS = np.linspace(-2, 2, 5)
POP_SIZE = 1e2


class Run:
    """Hold information relevant to a single Monte Carlo run."""

    def __init__(self, alpha, li, si):
        self.li = li
        self.si = si
        self.alpha = alpha
        self.survey_name = None
        self.ks_dm = None
        self.ks_snr = None
        self.ks_rate = None
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
        self.runs = None

    def set_up_surveys(self):
        """Set up surveys."""
        self.surveys = []
        for name in self.survey_names:
            survey = Survey(name=name)
            survey.set_beam(model='airy', n_sidelobes=1)
            self.surveys.append(survey)

    def generate(self, parallel=True):

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

                    pop.name = f'complex_alpha_{alpha}_lum_{li}_si_{si}'

                    for survey in self.surveys:
                        surv_pop = SurveyPopulation(pop, survey)
                        surv_pop.save()

                        # Save output
                        run = Run(alpha, li, si)
                        run.survey_name = survey.name
                        run.pop_name = surv_pop.name

                        # Add rate details
                        sr = surv_pop.source_rate
                        run.n_srcs = sr.det
                        run.n_days = sr.days
                        run.rate = sr.det / sr.days

                        # Add fit details
                        
                        runs.append(run)

            return runs

        if parallel:
            n_cpu = min([3, os.cpu_count() - 1])
            pprint(f'{os.cpu_count()} CPUs available')
            r = range(len(self.alphas))
            p = Parallel(n_jobs=n_cpu)(delayed(iter_alpha)(i) for i in tqdm(r))

        else:
            p = []
            for i in tqdm(range(len(self.alphas)), desc='Alphas'):
                p.append(iter_alpha(i))

        runs = [i for sublist in p for i in sublist]

    def get_runs(self):
        pass

    def plot_runs(self):
        pass


if __name__ == '__main__':
    mc = MonteCarlo(ALPHAS, LIS, SIS, SURVEY_NAMES, POP_SIZE)
    mc.generate()
