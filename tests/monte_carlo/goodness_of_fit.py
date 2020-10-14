from weighted_quantiles import median
from scipy.stats import ks_2samp
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt

from frbpoppy import unpickle, TNS, poisson_interval, pprint
from tests.rates.alpha_real import EXPECTED
from tests.convenience import plot_aa_style, rel_path
from simulations import SimulationOverview

NORM_SURV = 'parkes-htru'


class GoodnessOfFit:

    def __init__(self, runs=[1], calc=False, plot=False):
        self.runs = runs
        self.run_pars = {1: ['alpha', 'si', 'li'],
                         2: ['li', 'li_min', 'li_max'],
                         3: ['w_int_mean', 'w_int_std'],
                         4: ['dm_igm_slope', 'dm_host'],
                         5: ['alpha', 'si', 'li']}
        self.norm_surv = NORM_SURV
        self.so = SimulationOverview()
        self.tns = self.get_tns()

        for run in self.runs:
            if calc:
                self.calc(run)
            if plot:
                self.temp(run)
                self.plot(run)

    def get_tns(self):
        # Only get one-offs
        return TNS(repeaters=False, mute=True).df

    def dm(self, pop, survey_name):
        """Calculate GoodnessOfFit for DM distributions."""
        mask = (self.tns.survey == survey_name)
        try:
            ks_dm = ks_2samp(pop.frbs.dm, self.tns[mask].dm)[1]
        except ValueError:
            ks_dm = np.nan
        return ks_dm

    def snr(self, pop, survey_name):
        mask = (self.tns.survey == survey_name)
        try:
            ks_snr = ks_2samp(pop.frbs.snr, self.tns[mask].snr)[1]
        except ValueError:
            ks_snr = np.nan
        return ks_snr

    def rate(self, pop, survey_name, norm_uuid, run, errs=False):
        # Add rate details
        sr = pop.source_rate
        surv_sim_rate = sr.det / sr.days

        # Perhaps use at some stage
        if errs:
            p_int = poisson_interval(sr.det, sigma=1)
            surv_sim_rate_errs = [p/sr.days for p in p_int]

        # Don't bother if population is empty
        if sr.det == 0:
            return np.nan, 0

        # Determine ratio of detection rates
        if survey_name in EXPECTED:
            n_frbs, n_days = EXPECTED[survey_name]
        else:
            n_frbs, n_days = [np.nan, np.nan]
        surv_real_rate = n_frbs/n_days

        # Get normalisation properties
        norm_real_n_frbs, norm_real_n_days = EXPECTED[self.norm_surv]
        norm_pop = unpickle(f'mc/run_{run}/{norm_uuid}')
        norm_sim_n_frbs = norm_pop.source_rate.det
        norm_sim_n_days = norm_pop.source_rate.days
        norm_sim_rate = norm_sim_n_frbs / norm_sim_n_days
        norm_real_rate = norm_real_n_frbs / norm_real_n_days

        try:
            sim_ratio = surv_sim_rate / norm_sim_rate
            real_ratio = surv_real_rate / norm_real_rate
        except ZeroDivisionError:
            return np.nan, sr.det

        # TODO: Ideally update to include information on the Poisson intervals
        ks_rate = 1/np.abs(sim_ratio - real_ratio)
        if survey_name == self.norm_surv:
            ks_rate = np.nan

        return ks_rate, sr.det

    def calc(self, run):

        # For each requested run
        df = self.so.df

        pprint(f'Calculating goodness of fit for run {run}')

        # For each parameter
        for main_par in self.run_pars[run]:
            pprint(f' - {main_par}')
            pars = [e for e in self.run_pars[run] if e != main_par]
            pars.append('survey')

            # Loop over the remaining parameters
            for i, group in df[df.run == run].groupby(pars):
                s = group.survey.values[0]
                pprint(f'    - {list(zip(pars, i))}')

                # Calculate goodness of fit values for each simulation
                for row_ix, row in group.iterrows():
                    # Unpickle population
                    uuid = row.uuid
                    pop = unpickle(f'mc/run_{run}/{uuid}')
                    self.so.df.at[row_ix, 'dm_gof'] = self.dm(pop, s)
                    self.so.df.at[row_ix, 'snr_gof'] = self.snr(pop, s)

                    # Find corresponding rate normalisation population uuid
                    norm_mask = dict(zip(pars, i))
                    norm_mask['survey'] = self.norm_surv
                    norm_mask[main_par] = row[main_par]
                    k = norm_mask.keys()
                    v = norm_mask.values()
                    norm_uuid = df.loc[df[k].isin(v).all(axis=1), :].uuid
                    norm_uuid = norm_uuid.values[0]
                    rate_gof, n_det = self.rate(pop, s, norm_uuid, run)

                    # Get rate goodness of fit
                    self.so.df.at[row_ix, 'rate_gof'] = rate_gof
                    self.so.df.at[row_ix, 'n_det'] = n_det

        pprint(f'Saving the results for run {run}')

        # Calculate the weight per simluation
    def temp(self, run):
        df = self.so.df
        real = dict(self.tns.survey.value_counts())
        real_det = df.survey.map(real)
        self.so.df['weight'] = real_det + self.so.df['n_det']
        self.so.save()

    def plot(self, run):
        # Get data
        # For each requested run
        df = self.so.df

        # For each parameter
        for main_par in self.run_pars[run]:
            pprint(f'Plotting {main_par}')
            other_pars = [e for e in self.run_pars[run] if e != main_par]

            for compare_par in ['dm', 'snr', 'rate']:
                compare_col = f'{compare_par}_gof'

                pprint(f' - {compare_col}')
                for survey, group_surv in df[df.run == run].groupby('survey'):

                    pprint(f'    - {survey}')

                    # Set up plot
                    plot_aa_style()
                    plt.rcParams["figure.figsize"] = (5.75373*3, 5.75373*3)
                    plt.rcParams['figure.max_open_warning'] = 125
                    n_x = group_surv[other_pars[0]].nunique()
                    n_y = group_surv[other_pars[1]].nunique()
                    fig, ax = plt.subplots(n_x, n_y,
                                           sharex='col', sharey='row')
                    cmap = plt.get_cmap()

                    groups = group_surv.groupby(other_pars)
                    x = -1
                    for i, (other_pars_vals, group) in enumerate(groups):
                        bins = group[main_par].values
                        values = group[compare_col].values
                        bins, values = self.add_edges_to_hist(bins, values)
                        weight = np.median(group['weight'])

                        y = i % n_y
                        if y == 0:
                            x += 1
                        color = cmap(weight/group['weight'].max())
                        ax[y, x].step(bins, values, where='mid', color=color)
                        ax[y, x].set_title = str(other_pars_vals)

                        # Set axis label
                        if y == n_y - 1:
                            p = other_pars[0]
                            val = other_pars_vals[0]
                            ax[y, x].set_xlabel(f'{p} = {val:.2}')

                        if x == 0:
                            p = other_pars[1]
                            val = other_pars_vals[1]
                            ax[y, x].set_ylabel(f'{p} = {val:.2}')

                    for a in ax.flatten():
                        subset = df[df.run == run][main_par]
                        a.set_xlim(subset.min(), subset.max())

                    fig.suptitle(f'{main_par} - {compare_par} - {survey}')
                    plt.tight_layout()
                    plt.subplots_adjust(top=0.95)

                    # Save to subdirectory
                    path_to_save = rel_path(f'./plots/mc/{main_par}_run{run}/')
                    if not os.path.isdir(path_to_save):
                        os.mkdir(path_to_save)
                    path_to_save += f'{compare_par}_{survey}.pdf'
                    plt.savefig(path_to_save)
                    plt.clf()

    def add_edges_to_hist(self, bins, n, bin_type='lin'):
        """Add edges to histograms"""
        np.seterr(divide='ignore', invalid='ignore')

        if bin_type == 'lin':
            bin_dif = np.diff(bins)[-1]
            bins = np.insert(bins, 0, bins[0] - bin_dif)
            bins = np.insert(bins, len(bins), bins[-1] + bin_dif)
            n = np.insert(n, 0, 0)
            n = np.insert(n, len(n), 0)
        else:
            bin_dif = np.diff(np.log10(bins))[-1]
            bins = np.insert(bins, 0, 10**(np.log10(bins[0])-bin_dif))
            bins = np.insert(bins, len(bins), 10**(np.log10(bins[-1])+bin_dif))
            n = np.insert(n, 0, 0)
            n = np.insert(n, len(n), 0)
        return bins, n

    def weighted_median(self, df):
        gofs = pd.melt(df[['dm_gof', 'snr_gof', 'rate_gof']]).value.values
        weights = np.concatenate(df['weight']*3)
        return median(gofs, weights)


if __name__ == '__main__':
    GoodnessOfFit(runs=[1], calc=False, plot=True)
