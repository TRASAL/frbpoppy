from weighted_quantiles import median
from scipy.stats import ks_2samp
import numpy as np
import os
import matplotlib.pyplot as plt

from frbpoppy import unpickle, TNS, poisson_interval, pprint
from tests.rates.alpha_real import EXPECTED
from tests.convenience import plot_aa_style, rel_path
from simulations import SimulationOverview, POP_SIZE

NORM_SURV = 'parkes-htru'


class GoodnessOfFit:

    def __init__(self):
        self.run_pars = {1: ['alpha', 'si', 'li'],
                         2: ['li', 'lum_min', 'lum_max'],
                         3: ['w_mean', 'w_std'],
                         4: ['dm_igm_slope', 'dm_host']}
        self.norm_surv = NORM_SURV
        self.so = SimulationOverview()
        self.tns = self.get_tns()

    def get_tns(self):
        # Only get one-offs
        return TNS(repeaters=False, mute=True, update=False).df

    def dm(self, pop, survey_name):
        """Calculate GoodnessOfFit for DM distributions."""
        mask = ((self.tns.survey == survey_name) & (self.tns.dm <= 950))
        try:
            ks_dm = ks_2samp(pop.frbs.dm, self.tns[mask].dm)[1]
        except ValueError:
            ks_dm = np.nan
        return ks_dm

    def snr(self, pop, survey_name):
        mask = ((self.tns.survey == survey_name) & (self.tns.dm <= 950))
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

        if norm_sim_rate == 0:
            norm_sim_rate = POP_SIZE / norm_sim_n_days

        sim_ratio = surv_sim_rate / norm_sim_rate
        real_ratio = surv_real_rate / norm_real_rate

        diff = np.abs(sim_ratio - real_ratio)
        if diff == 0:
            rate_diff = 1e-3
        else:
            rate_diff = 1 / diff

        return rate_diff, pop.n_sources()

    def calc_gofs(self, run):

        # For each requested run
        self.so = SimulationOverview()
        par_set = self.so.df[self.so.df.run == run].par_set.iloc[0]
        pprint(f'Calculating goodness of fit for run {run}, par set {par_set}')
        pars = self.run_pars[par_set]
        values = []

        # Loop through all combination of parameters
        for values, group in self.so.df[self.so.df.run == run].groupby(pars):
            pprint(f'    - {list(zip(pars, values))}')
            # Calculate goodness of fit values for each simulation
            for row_ix, row in group.iterrows():
                survey_name = row.survey
                uuid = row.uuid
                pop = unpickle(f'mc/run_{run}/{uuid}')

                # Apply a DM cutoff
                mask = (pop.frbs.dm <= 950)
                pop.frbs.apply(mask)
                pop.source_rate.det = pop.n_sources() * pop.source_rate.f_area

                dm_gof = self.dm(pop, survey_name)
                snr_gof = self.snr(pop, survey_name)
                self.so.df.at[row_ix, 'dm_gof'] = dm_gof
                self.so.df.at[row_ix, 'snr_gof'] = snr_gof

                if pop.n_sources() == 0:
                    self.so.df.at[row_ix, 'weight'] = 0
                    self.so.df.at[row_ix, 'n_det'] = pop.n_sources()
                    pprint(f'        -  No sources in {survey_name}')
                    continue

                # Find corresponding rate normalisation population uuid
                norm_mask = dict(zip(pars, values))
                norm_mask['survey'] = self.norm_surv
                norm_mask['run'] = run
                k = norm_mask.keys()
                v = norm_mask.values()
                norm_uuid = group.loc[group[k].isin(v).all(axis=1), :].uuid
                norm_uuid = norm_uuid.values[0]
                rate_diff, n_det = self.rate(pop, survey_name, norm_uuid, run)

                # Get rate weighting
                self.so.df.at[row_ix, 'weight'] = rate_diff
                self.so.df.at[row_ix, 'n_det'] = n_det

        pprint(f'Saving the results for run {run}')
        # Best matching in terms of rates
        max_w = np.nanmax(self.so.df.weight)
        self.so.df.loc[self.so.df.weight == 1e3]['weight'] = max_w
        self.so.save()

    def plot(self, run):
        # Get data
        # For each requested run
        df = self.so.df
        par_set = df[df.run == run].par_set.iloc[0]

        # For each parameter
        for main_par in self.run_pars[par_set]:
            pprint(f'Plotting {main_par}')
            other_pars = [e for e in self.run_pars[par_set] if e != main_par]

            for compare_par in ['dm', 'snr']:
                compare_col = f'{compare_par}_gof'

                pprint(f' - {compare_col}')
                for survey, group_surv in df[df.run == run].groupby('survey'):

                    pprint(f'    - {survey}')

                    # Set up plot
                    plot_aa_style()
                    plt.rcParams["figure.figsize"] = (5.75373*3, 5.75373*3)
                    plt.rcParams['figure.max_open_warning'] = 125
                    n_x = group_surv[other_pars[0]].nunique()
                    if len(other_pars) > 1:
                        n_y = group_surv[other_pars[1]].nunique()
                    else:
                        n_y = 1
                    fig, ax = plt.subplots(n_x, n_y,
                                           sharex='col', sharey='row')

                    groups = group_surv.groupby(other_pars)
                    x = -1
                    for i, (other_pars_vals, group) in enumerate(groups):
                        bins = group[main_par].values
                        values = group[compare_col].values
                        bins, values = self.add_edges_to_hist(bins, values)

                        if n_y > 1:
                            y = i % n_y
                            if y == 0:
                                x += 1
                            a = ax[y, x]
                        else:
                            y = i
                            a = ax[y]

                        a.step(bins, values, where='mid')
                        a.set_title = str(other_pars_vals)

                        diff = np.diff(bins)
                        if diff[1] != diff[0]:
                            a.set_xscale('log')

                        # Set axis label
                        if y == n_y - 1:
                            p = other_pars[0]
                            if isinstance(other_pars_vals, float):
                                val = other_pars_vals
                            else:
                                val = other_pars_vals[0]
                            p = p.replace('_', ' ')
                            a.set_xlabel(f'{p} = {val:.2}')

                        if x == 0:
                            p = other_pars[1]
                            val = other_pars_vals[1]
                            p = p.replace('_', ' ')
                            a.set_ylabel(f'{p} = {val:.2}')

                    # Set axis limits
                    subset = df[df.run == run][main_par]
                    y_subset = group_surv[compare_col].copy()
                    try:
                        low = np.nanmin(y_subset)
                        high = np.nanmax(y_subset)
                    except ValueError:
                        low = 0.0001
                        high = 1
                    log = False
                    if low > 0 and high > 0:
                        log = True

                    for a in ax.flatten():
                        a.set_xlim(subset.min(), subset.max())
                        if log:
                            a.set_yscale('log', nonposy='clip')
                            a.set_ylim(low, high)

                    p = main_par.replace('_', ' ')
                    fig.suptitle(f'{p} - {compare_par} - {survey}')
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
        else:
            bin_dif = np.diff(np.log10(bins))[-1]
            bins = np.insert(bins, 0, 10**(np.log10(bins[0])-bin_dif))
            bins = np.insert(bins, len(bins), 10**(np.log10(bins[-1])+bin_dif))

        n = np.insert(n, 0, np.nan)
        n = np.insert(n, len(n), np.nan)

        return bins, n

    def weighted_median(self, df):
        dm_gof = df['dm_gof'].values
        dm_weight = df['n_det'].values
        snr_gof = df['snr_gof'].values
        snr_weight = df['n_det'].values
        gofs = np.concatenate([dm_gof, snr_gof])
        weights = np.concatenate([dm_weight, snr_weight])
        return median(gofs, weights)

    def calc_global_max(self, run):
        self.so = SimulationOverview()
        df = self.so.df[self.so.df.run == run]
        par_set = df[df.run == run].par_set.iloc[0]
        cols = self.run_pars[par_set]
        values = []
        gofs = []

        # Loop through all combination of parameters
        for value, group in df.groupby(cols):
            gof = self.weighted_median(group)
            values.append(value)
            gofs.append(gof)

        gofs = np.array(gofs)

        # Find maximum (best gof)
        if np.isnan(gofs).all():
            return dict(zip(cols, [(np.nan, np.nan) for i in cols]))
        else:
            best_ix = np.nanargmax(gofs)
            best_values = values[best_ix]
            best_gofs = [gofs[best_ix]]*len(cols)

        return dict(zip(cols, zip(best_values, best_gofs)))
