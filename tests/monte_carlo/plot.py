from frbpoppy import hist
from goodness_of_fit import GoodnessOfFit, Fit
from matplotlib.lines import Line2D
from tests.convenience import plot_aa_style, rel_path
import matplotlib.pyplot as plt
import numpy as np


class Plot():
    """Plot runs."""

    def __init__(self):
        plot_aa_style()
        plt.rcParams['figure.figsize'] = (5.75373, (5.75373/3)*4)
        plt.rcParams['font.size'] = 9
        plt.rcParams['xtick.major.pad'] = 10
        plt.rcParams['ytick.major.pad'] = 10
        plt.rcParams['axes.titlepad'] = 10
        self.fig, self.axes = plt.subplots(4, 3, sharey='row')
        self.colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

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

        # Testing
        self.gen_test()

        # Plot run rectangles
        self.runs()

        # Save plot
        plt.tight_layout(rect=[0, 0, 0.98, 1])
        plt.savefig(rel_path('./plots/mc.pdf'))

    def alpha(self):
        ax = self.axes[0, 0]
        ax.set_yscale('log')
        ax.set_ylabel(r'g.o.f.')
        parm = r'\alpha'
        ax.set_xlabel(rf'${parm}$')

        best_value, best_value_err = 1, 1
        ax.set_title(rf'${parm}={best_value}\pm{best_value_err}$',
                     color=self.colors[2])

    def si(self):
        ax = self.axes[0, 1]
        parm = r'\text{si}'
        ax.set_xlabel(rf'${parm}$')
        best_value, best_value_err = 1, 1
        ax.set_title(rf'${parm}={best_value}\pm{best_value_err}$',
                     color=self.colors[2])

    def li(self):
        ax = self.axes[0, 2]
        parm = r'\text{lum$_{\text{i}}$}'
        ax.set_xlabel(rf'${parm}$')
        best_value, best_value_err = 1, 1
        ax.set_title(rf'${parm}={best_value}\pm{best_value_err}$',
                     color=self.colors[2])

    def li_2(self):
        ax = self.axes[1, 0]
        ax.set_yscale('log')
        ax.set_ylabel(r'g.o.f.')
        parm = r'\text{lum$_{\text{i}}$}'
        ax.set_xlabel(rf'${parm}$')
        best_value, best_value_err = 1, 1
        ax.set_title(rf'${parm}={best_value}\pm{best_value_err}$',
                     color=self.colors[2])

    def lum_min(self):
        ax = self.axes[1, 1]
        parm = r'\text{lum$_{\text{min}}$}'
        ax.set_xlabel(rf'${parm}$')
        best_value, best_value_err = 1, 1
        ax.set_title(rf'${parm}={best_value}\pm{best_value_err}$',
                     color=self.colors[2])

    def lum_max(self):
        ax = self.axes[1, 2]
        parm = r'\text{lum$_{\text{max}}$}'
        ax.set_xlabel(rf'${parm}$')
        best_value, best_value_err = 1, 1
        ax.set_title(rf'${parm}={best_value}\pm{best_value_err}$',
                     color=self.colors[2])

    def w_int_mean(self):
        ax = self.axes[2, 0]
        ax.set_yscale('log')
        ax.set_ylabel(r'g.o.f.')
        parm = r'\text{w$_{\text{int, }\mu}$}'
        ax.set_xlabel(rf'${parm}$')
        best_value, best_value_err = 1, 1
        ax.set_title(rf'${parm}={best_value}\pm{best_value_err}$',
                     color=self.colors[2])

    def w_int_std(self):
        ax = self.axes[2, 1]
        parm = r'\text{w$_{\text{int, }\sigma}$}'
        ax.set_xlabel(rf'${parm}$')
        best_value, best_value_err = 1, 1
        ax.set_title(rf'${parm}={best_value}\pm{best_value_err}$',
                     color=self.colors[2])

    def dm_igm_slope(self):
        ax = self.axes[3, 0]
        ax.set_yscale('log')
        ax.set_ylabel(r'g.o.f.')
        parm = r'\text{DM$_{\text{IGM, slope}}$}'
        ax.set_xlabel(rf'${parm}$')
        best_value, best_value_err = 1, 1
        ax.set_title(rf'${parm}={best_value}\pm{best_value_err}$',
                     color=self.colors[2])

    def dm_host(self):
        ax = self.axes[3, 1]
        parm = r'\text{DM$_{\text{Host}}$}'
        ax.set_xlabel(rf'${parm}$')
        best_value, best_value_err = 1, 1
        ax.set_title(rf'${parm}={best_value}\pm{best_value_err}$',
                     color=self.colors[2])

    def legend(self):
        ax = self.axes[2, 2]

        # Add legend elements
        elements = []

        line = Line2D([0], [0], color=self.colors[0])
        label = 'Simulated'
        elements.append((line, label))

        line = Line2D([0], [0], color=self.colors[1])
        label = 'Check'
        elements.append((line, label))

        line = Line2D([0], [0], color=self.colors[3])
        label = 'Preset'
        elements.append((line, label))

        line = Line2D([0], [0], color=self.colors[2])
        label = 'Best fit'
        elements.append((line, label))

        lines, labels = zip(*elements)
        self.fig.legend(lines, labels, bbox_to_anchor=(0.80, 0.4),
                        loc='center')

        ax.set_axis_off()

    def runs(self):
        sw = 0.280  # Subplot width
        sh = 0.245  # Subplot height
        vp = 0.004  # Vertical padding
        le = 0.655  # Left coordinate
        b = 0.8075  # Bottom coordinate
        self.make_rectangle(le-2*sw, b-vp, 0.305*3-0.03, 0.147+2*vp,
                            'Run 1, 5')
        self.make_rectangle(le-2*sw, b-vp-sh, 0.305*3-0.03, 0.147+2*vp,
                            'Run 2')

        sh += 0.001
        self.make_rectangle(le-2*sw, b-vp-2*sh, 0.305*2-0.02, 0.147+2*vp,
                            'Run 3')
        self.make_rectangle(le-sw-sw, b-vp-3*sh, 0.305*2-0.02, 0.147+2*vp,
                            'Run 4')

    def make_rectangle(self, left, bottom, width, height, name, alpha=1):
        rect = plt.Rectangle((left, bottom), width, height, fill=False,
                             color="k", lw=0.5, zorder=900, alpha=alpha,
                             transform=self.fig.transFigure, figure=self.fig)
        rect.set_transform(self.fig.transFigure)
        right = left + width
        top = bottom + height
        self.fig.text(right, 0.5 * (bottom + top), name,
                      horizontalalignment='center',
                      verticalalignment='center',
                      rotation='vertical',
                      transform=self.fig.transFigure,
                      zorder=1000,
                      backgroundcolor='white',
                      color='k',
                      alpha=alpha,
                      fontsize=7)
        self.fig.patches.extend([rect])

    def gen_test(self):
        for i, ax in enumerate(self.axes.flatten()):
            if i == 8 or i == 11:
                continue
            # simluation data
            data = np.random.normal(0, 0.1, 1000)
            bin_centres, values = hist(data, bins=11)
            ax.step(bin_centres, 10**values, where='mid')

            # convergance check
            data = np.random.normal(0, 0.05, 1000)
            bin_centres, values = hist(data, bins=11)
            if i < 3:
                ax.step(bin_centres, 10**values, where='mid')

            # Preset value
            ax.axvline(x=0.1, color=self.colors[3], lw=1)

            # Postset value
            ax.axvline(x=0, color=self.colors[2], lw=1)


if __name__ == '__main__':
    Plot()
