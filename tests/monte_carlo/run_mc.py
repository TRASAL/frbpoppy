"""Link together all classes to run a full Monte Carlo."""
import pandas as pd
import numpy as np
import frbpoppy.paths
import os

from simulations import MonteCarlo, POP_SIZE
from goodness_of_fit import GoodnessOfFit
from plot import Plot

GENERATE = True
CALC_GOFS = True
RUNS = [10]


class RunOverview:
    """Gather input for each run."""

    def __init__(self, load_csv=True):
        p = frbpoppy.paths.populations()
        self.filename = f'{p}mc/run_overview.csv'

        if load_csv and os.path.isfile(self.filename):
            self.df = self.load()
        else:
            self.df = self.gen_runs()

    def gen_run(self):
        return {'alpha': None,
                'si': None,
                'li': None,
                'lum_min': None,
                'lum_max': None,
                'w_mean': None,
                'w_std': None,
                'dm_igm_slope': None,
                'dm_host': None,
                'execute': True,
                'par_set': 0,
                'run': 0}

    def gen_runs(self):
        runs = []
        for i in range(10):
            r = self.gen_run()
            r['run_number'] = i + 1
            r['execute'] = True
            r['par_set'] = i % 4 + 1
            if i == 9:  # Holder for best values
                r['execute'] = False
            runs.append(r)
        df = pd.DataFrame(runs)
        df.set_index('run', inplace=True)
        return df

    def load(self):
        df = pd.read_csv(self.filename)
        df.run = df.run.astype(int)
        df.par_set = df.par_set.astype(int)
        df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
        return df

    def save(self, df=None):
        if df is None:
            df = self.df
        df.to_csv(self.filename)


if __name__ == '__main__':
    print('Commencing')
    runs = RunOverview(load_csv=True)
    mc = MonteCarlo(pop_size=POP_SIZE)

    for i, run in runs.df.iterrows():

        run = runs.df.iloc[i]

        print('='*50)
        print(f'On Run {run.run} with par_set {run.par_set}')
        print('='*50)
        print(run)

        if run.run not in RUNS:
            continue

        # Generate parameter sets
        if GENERATE:
            if run.par_set == 1:
                mc.gen_par_set_1(lum_min=run.lum_min,
                                 lum_max=run.lum_max,
                                 w_mean=run.w_mean,
                                 w_std=run.w_std,
                                 dm_igm_slope=run.dm_igm_slope,
                                 dm_host=run.dm_host,
                                 run=run.run)
            if run.par_set == 2:
                mc.gen_par_set_2(alpha=run.alpha,
                                 si=run.si,
                                 w_mean=run.w_mean,
                                 w_std=run.w_std,
                                 dm_igm_slope=run.dm_igm_slope,
                                 dm_host=run.dm_host,
                                 run=run.run)
            if run.par_set == 3:
                mc.gen_par_set_3(alpha=run.alpha,
                                 si=run.si,
                                 li=run.li,
                                 lum_min=run.lum_min,
                                 lum_max=run.lum_max,
                                 dm_igm_slope=run.dm_igm_slope,
                                 dm_host=run.dm_host,
                                 run=run.run)
            if run.par_set == 4:
                mc.gen_par_set_4(alpha=run.alpha,
                                 si=run.si,
                                 li=run.li,
                                 lum_min=run.lum_min,
                                 lum_max=run.lum_max,
                                 w_mean=run.w_mean,
                                 w_std=run.w_std,
                                 run=run.run)

        # Determine the goodness of fit
        gf = GoodnessOfFit()
        if CALC_GOFS:
            gf.calc_gofs(run.run)

        # Find global maximums
        gf = GoodnessOfFit()
        gms = gf.calc_global_max(run.run)
        print('\n')
        print(f'   Best fits from run {run.run}-> {gms}')
        print('\n')

        # Adapt the input for future runs
        for j in range(i+1, len(runs.df)):
            for par in gms:
                runs.df.at[j, par] = gms[par][0]

        runs.save()

        Plot()
