from frbpoppy import hist
from goodness_of_fit import GoodnessOfFit
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit
from tests.convenience import plot_aa_style, rel_path
import matplotlib.pyplot as plt
import numpy as np


class Plot():
    """Plot runs."""

    def __init__(self):
        plot_aa_style()
        plt.rcParams['figure.figsize'] = (5.75373, (5.75373/3)*4)
        plt.rcParams['font.size'] = 9
        # plt.rcParams['xtick.major.pad'] = 10
        # plt.rcParams['ytick.major.pad'] = 10
        # plt.rcParams['axes.titlepad'] = 10
        self.fig, self.axes = plt.subplots(4, 3, sharey='row')
        self.colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        self.gf = GoodnessOfFit()
        self.df = self.gf.so.df

        # Calculate global maximums
        self.gm = {}
        for run in self.df.run.unique():
            self.gm[run] = self.gf.calc_global_max(run)
        print(self.gm)

        # Plot various subplots
        self.alpha()
        self.si()
        self.li()
        self.li_2()
        self.lum_min()
        self.lum_max()
        self.w_int_mean()
        self.w_int_std()
        self.legend()
        self.dm_igm_slope()
        self.dm_host()
        self.axes[3, 2].set_axis_off()

        # Plot run rectangles
        # self.runs()

        # Save plot
        plt.tight_layout()  # rect=[0, 0, 0.98, 1])
        plt.subplots_adjust(wspace=0.1)
        plt.savefig(rel_path('./plots/mc/mc.pdf'))

    def alpha(self):
        ax = self.axes[0, 0]
        parm = 'alpha'
        label = r'$\alpha$'
        runs = [1, 5, 8]

        ax.set_yscale('log', nonposy='clip')
        ax.set_ylabel(r'GoF')
        ax.set_xlabel(label)
        best_value = np.nan

        # Plot runs
        for i, run in enumerate(runs):
            bins, gofs = self.get_data(run, parm)
            ax.step(bins, gofs, where='mid')

            # Plot global maximum
            try:
                best_value, best_gof = self.gm[run][parm]
                ax.plot([best_value]*2, [1e-1, best_gof], color=self.colors[i], linestyle='--')
                ax.scatter([best_value], [best_gof], marker='x', color=self.colors[i])
            except KeyError:
                i -= 1
                continue

            if i==1 and not np.isnan(best_value):
                title = fr'{label}=${best_value:.1f}$'
                ax.set_title(title, fontsize=10, color=self.colors[i])

    def si(self):
        ax = self.axes[0, 1]
        parm = 'si'
        label = r'\text{si}'
        runs = [1, 5, 8]

        ax.set_yscale('log', nonposy='clip')
        ax.set_xlabel(label)
        best_value = np.nan

        # Plot runs
        for i, run in enumerate(runs):
            bins, gofs = self.get_data(run, parm)
            ax.step(bins, gofs, where='mid')

            # Plot global maximum
            try:
                best_value, best_gof = self.gm[run][parm]
                ax.plot([best_value]*2, [1e-1, best_gof], color=self.colors[i], linestyle='--')
                ax.scatter([best_value], [best_gof], marker='x', color=self.colors[i])
            except KeyError:
                i -= 1
                continue

            if i == 1 and not np.isnan(best_value):
                title = fr'{label}=${best_value:.1f}$'
                ax.set_title(title, fontsize=10, color=self.colors[i])

    def li(self):
        ax = self.axes[0, 2]
        parm = 'li'
        # label = r'\text{lum$_{\text{i}}$}'
        label = 'li'
        runs = [1, 5, 8]

        ax.set_yscale('log', nonposy='clip')
        ax.set_xlabel(label)
        ax_right = ax.twinx()
        ax_right.set_ylabel('Set 1', labelpad=10)
        ax_right.tick_params(axis='y', which='both', right=False, labelright=False)
        best_value = np.nan

        # Plot runs
        for i, run in enumerate(runs):
            bins, gofs = self.get_data(run, parm)
            ax.step(bins, gofs, where='mid')

            # Plot global maximum
            try:
                best_value, best_gof = self.gm[run][parm]
                ax.plot([best_value]*2, [1e-1, best_gof], color=self.colors[i], linestyle='--')
                ax.scatter([best_value], [best_gof], marker='x', color=self.colors[i])
            except KeyError:
                i -= 1
                continue

            if i == 1 and not np.isnan(best_value):
                title = fr'{label}=${best_value:.1f}$'
                ax.set_title(title, fontsize=10, color=self.colors[i])

    def li_2(self):
        ax = self.axes[1, 0]
        ax.set_ylabel(r'GoF')
        parm = 'li'
        label = 'li'
        runs = [2]

        ax.set_yscale('log', nonposy='clip')
        ax.set_xlabel(label)
        best_value = np.nan

        # Plot runs
        for i, run in enumerate(runs):
            bins, gofs = self.get_data(run, parm)
            ax.step(bins, gofs, where='mid')

            # Plot global maximum
            try:
                best_value, best_gof = self.gm[run][parm]
                ax.plot([best_value]*2, [1e-1, best_gof], color=self.colors[i], linestyle='--')
                ax.scatter([best_value], [best_gof], marker='x', color=self.colors[i])
            except KeyError:
                i -= 1
                continue

        # if not np.isnan(best_value):
        #     title = fr'{label}=${best_value:.1f}$'
        #     ax.set_title(title, fontsize=10, color=self.colors[i])

    def lum_min(self):
        ax = self.axes[1, 1]
        ax.set_xscale('log')
        label = r'\text{lum$_{\text{min}}$}'
        parm = 'lum_min'
        runs = [2]

        ax.set_yscale('log', nonposy='clip')
        ax.set_xlabel(label + r' (erg s$^{-1}$)')
        best_value = np.nan

        # Plot runs
        for i, run in enumerate(runs):
            bins, gofs = self.get_data(run, parm)
            ax.step(bins, gofs, where='mid')

            # Plot global maximum
            try:
                best_value, best_gof = self.gm[run][parm]
                ax.plot([best_value]*2, [1e-1, best_gof], color=self.colors[i], linestyle='--')
                ax.scatter([best_value], [best_gof], marker='x', color=self.colors[i])
            except KeyError:
                i -= 1
                continue

        # if not np.isnan(best_value):
        #     title = fr'{label}=${best_value:.1e}$'
        #     ax.set_title(title, fontsize=10, color=self.colors[i])

    def lum_max(self):
        ax = self.axes[1, 2]
        ax.set_xscale('log')
        label = r'\text{lum$_{\text{max}}$}'
        parm = 'lum_max'
        runs = [2]

        ax.set_yscale('log', nonposy='clip')
        ax.set_xlabel(label + r' (erg s$^{-1}$)')
        ax_right = ax.twinx()
        ax_right.set_ylabel('Set 2', labelpad=10)
        ax_right.tick_params(axis='y', which='both', right=False, labelright=False)
        best_value = np.nan

        # Plot runs
        for i, run in enumerate(runs):
            bins, gofs = self.get_data(run, parm)
            ax.step(bins, gofs, where='mid')

            # Plot global maximum
            try:
                best_value, best_gof = self.gm[run][parm]
                ax.plot([best_value]*2, [1e-1, best_gof], color=self.colors[i], linestyle='--')
                ax.scatter([best_value], [best_gof], marker='x', color=self.colors[i])
            except KeyError:
                i -= 1
                continue

        # if not np.isnan(best_value):
        #     title = fr'{label}=${best_value:.1e}$'
        #     ax.set_title(title, fontsize=10, color=self.colors[i])

    def w_int_mean(self):
        ax = self.axes[2, 0]
        ax.set_xscale('log')
        ax.set_ylabel(r'GoF')
        label = r'\text{w$_{\text{int, mean}}$}'
        parm = 'w_mean'
        runs = [3, 6, 9]

        ax.set_yscale('log', nonposy='clip')
        ax.set_xlabel(fr'{label} (ms)')
        best_value = np.nan

        # Plot runs
        for i, run in enumerate(runs):
            bins, gofs = self.get_data(run, parm)
            ax.step(bins, gofs, where='mid')

            # Plot global maximum
            try:
                best_value, best_gof = self.gm[run][parm]
                ax.plot([best_value]*2, [1e-1, best_gof], color=self.colors[i], linestyle='--')
                ax.scatter([best_value], [best_gof], marker='x', color=self.colors[i])
            except KeyError:
                i -= 1
                continue

            if i == 1 and not np.isnan(best_value):
                title = fr'{label}=${best_value:.1e}$'
                ax.set_title(title, fontsize=10, color=self.colors[i])

    def w_int_std(self):
        ax = self.axes[2, 1]
        label = r'\text{w$_{\text{int, std}}$}'
        parm = 'w_std'
        runs = [3, 6, 9]

        ax.set_yscale('log', nonposy='clip')
        ax.set_xlabel(fr'{label} (ms)')
        ax_right = ax.twinx()
        ax_right.set_ylabel('Set 3', labelpad=10)
        ax_right.tick_params(axis='y', which='both', right=False, labelright=False)
        best_value = np.nan

        # Plot runs
        for i, run in enumerate(runs):
            bins, gofs = self.get_data(run, parm)
            ax.step(bins, gofs, where='mid')

            # Plot global maximum
            try:
                best_value, best_gof = self.gm[run][parm]
                ax.plot([best_value]*2, [1e-1, best_gof], color=self.colors[i], linestyle='--')
                ax.scatter([best_value], [best_gof], marker='x', color=self.colors[i])
            except KeyError:
                i -= 1
                continue

            if i == 1 and not np.isnan(best_value):
                title = fr'{label}=${best_value:.1f}$'
                ax.set_title(title, fontsize=10, color=self.colors[i])

    def dm_igm_slope(self):
        ax = self.axes[3, 0]
        ax.set_ylabel(r'GoF')
        label = r'\text{DM$_{\text{IGM, slope}}$}'
        parm = 'dm_igm_slope'
        runs = [4, 7, 10]

        ax.set_yscale('log', nonposy='clip')
        ax.set_xlabel(label + r' ($\textrm{pc}\ \textrm{cm}^{-3}$)')
        best_value = np.nan

        # Plot runs
        for i, run in enumerate(runs):
            bins, gofs = self.get_data(run, parm)
            ax.step(bins, gofs, where='mid')

            # Plot global maximum
            try:
                best_value, best_gof = self.gm[run][parm]
                ax.plot([best_value]*2, [1e-1, best_gof], color=self.colors[i], linestyle='--')
                ax.scatter([best_value], [best_gof], marker='x', color=self.colors[i])
            except KeyError:
                i -= 1
                continue

            if i == 1 and not np.isnan(best_value):
                title = fr'{label}=${best_value:.0f}$'
                ax.set_title(title, fontsize=10, color=self.colors[i])

    def dm_host(self):
        ax = self.axes[3, 1]
        label = r'\text{DM$_{\text{Host}}$}'
        parm = 'dm_host'
        runs = [4, 7, 10]

        ax.set_yscale('log', nonposy='clip')
        ax.set_xlabel(label + r' ($\textrm{pc}\ \textrm{cm}^{-3}$)')
        ax_right = ax.twinx()
        ax_right.set_ylabel('Set 4', labelpad=10)
        ax_right.tick_params(axis='y', which='both', right=False, labelright=False)
        best_value = np.nan

        # Plot runs
        for i, run in enumerate(runs):
            bins, gofs = self.get_data(run, parm)
            ax.step(bins, gofs, where='mid')

            # Plot global maximum
            try:
                best_value, best_gof = self.gm[run][parm]
                ax.plot([best_value]*2, [1e-1, best_gof], color=self.colors[i], linestyle='--')
                ax.scatter([best_value], [best_gof], marker='x', color=self.colors[i])
            except KeyError:
                i -= 1
                continue

            if i == 1 and not np.isnan(best_value):
                title = fr'{label}=${best_value:.0f}$'
                ax.set_title(title, fontsize=10, color=self.colors[i])

    def legend(self):
        ax = self.axes[2, 2]

        # Add legend elements
        elements = []

        line = Line2D([0], [0], color=self.colors[0])
        elements.append((line, r'Cycle 1'))

        line = Line2D([0], [0], color=self.colors[1])
        elements.append((line, r'Cycle 2'))

        line = Line2D([0], [0], color=self.colors[2])
        elements.append((line, r'Cycle 3'))

        line = Line2D([0], [0], color='grey', linestyle='--')
        label = r'Max GoF'
        elements.append((line, label))

        lines, labels = zip(*elements)
        self.fig.legend(lines, labels, bbox_to_anchor=(0.84, 0.4),
                        loc='center')

        ax.set_axis_off()

    def get_data(self, run_number, par):
        df = self.df[self.df.run == run_number]
        if df.empty:
            return [np.nan], [np.nan]
        gofs = []
        bins = []
        for bin_val, group in df.groupby(par):
            gof = self.gf.weighted_median(group)
            gofs.append(gof)
            bins.append(bin_val)
        bins = np.array(bins)
        gofs = np.array(gofs)
        diff = np.diff(bins)
        bin_type = 'lin'
        if not np.isclose(diff[0], diff[1]):
            bin_type = 'log'
        bins, gofs = self.gf.add_edges_to_hist(bins, gofs, bin_type=bin_type)
        gofs[np.isnan(gofs)] = 1e-1
        return bins, gofs


if __name__ == '__main__':
    Plot()
