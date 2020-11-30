"""Check the log N log F slope for future surveys."""
import numpy as np
import matplotlib.pyplot as plt
from copy import copy

from frbpoppy import CosmicPopulation, Survey, LargePopulation, SurveyPopulation, hist
from frbpoppy import unpickle, pprint
import frbpoppy.direction_dists as did
import frbpoppy.galacticops as go

from tests.convenience import plot_aa_style, rel_path

from tests.rates.alpha_real import EXPECTED

MAKE = True
SURVEYS = ('parkes-htru',
           'wsrt-apertif',
           'fast-crafts',
           'puma-full',
           'chord',
           'ska1-low',
           'ska1-mid')
SIZE = 5e4


if MAKE:
    # Calculate the fraction of the sky that the survey covers
    surv_f_area = {}
    for name in SURVEYS:
        pop = CosmicPopulation.simple(5e5)
        pop.gen_direction()
        survey = Survey(name)
        mask = survey.in_region(pop.frbs.ra, pop.frbs.dec,
                                pop.frbs.gl, pop.frbs.gb)
        in_surv_region = np.sum(mask)
        tot_region = len(mask)

        area_sky = 4*np.pi*(180/np.pi)**2   # In sq. degrees
        f_area = (survey.beam_size/area_sky)*(tot_region/in_surv_region)

        surv_f_area[name] = f_area
        print(f'{name} covers {f_area*100}% of the sky')

    surv_pops = []

    for name in SURVEYS:
        # Set up survey
        survey = Survey(name)
        if name in ('parkes-htru', 'wsrt-apertif'):
            survey.set_beam(model=name)

        # Set up CosmicPopulation
        pop = CosmicPopulation.optimal(SIZE, generate=False)

        # Only generate FRBs in the survey region
        pop.set_direction(model='uniform',
                          min_ra=survey.ra_min,
                          max_ra=survey.ra_max,
                          min_dec=survey.dec_min,
                          max_dec=survey.dec_max)

        # Parkes also has galactic limits:
        if name == 'parkes-htru':
            pop.gen_index()
            pop.gen_dist()
            pop.gen_time()

            # Generate FRBs just within the galactic constraints
            pop.gen_direction()

            # Gather ra, dec coordinate limits
            lims = {'min_ra': survey.ra_min, 'max_ra': survey.ra_max,
                    'min_dec': survey.dec_min, 'max_dec': survey.dec_max}

            def sample(n_gen):
                ra, dec = did.uniform(n_srcs=n_gen, **lims)
                gl, gb = go.radec_to_lb(ra, dec, frac=True)
                coords = [ra, dec, gl, gb]
                return coords

            def accept(coords):
                return survey.in_region(*coords)

            coords = sample(int(SIZE))
            mask = accept(coords)
            reject, = np.where(~mask)
            while reject.size > 0:
                fill = sample(reject.size)
                mask = accept(fill)
                for i in range(len(coords)):
                    coords[i][reject[mask]] = fill[i][mask]
                reject = reject[~mask]

            # Assign the values
            frbs = pop.frbs
            frbs.ra, frbs.dec = coords[0], coords[1]
            frbs.gl, frbs.gb = coords[2], coords[3]

            # Continue with generation
            pop.gen_gal_coords()
            pop.gen_dm()
            pop.gen_w()
            pop.gen_lum()
            pop.gen_si()
        else:
            pop.generate()

        surv_pop = SurveyPopulation(pop, survey, scale_by_area=False)
        surv_pop.source_rate.f_area = surv_f_area[name]
        surv_pop.source_rate.scale_by_area()
        # surv_pop.save()
        surv_pops.append(surv_pop)

else:
    surv_pops = []
    for name in SURVEYS:
        surv_pops.append(unpickle(f'optimal_{name}'))

# Start plot
plot_aa_style(cols=2)
plt.rcParams["figure.figsize"] = (3.556*3, 3.556)
fig, axes = plt.subplots(1, 3)

for ax in axes.flatten():
    ax.set_aspect('auto')

# Get norm pop
y = 0
ys = []
names = []
rates = []
norm_sim_rate = surv_pops[0].source_rate.det
norm_real_rate = EXPECTED['parkes-htru'][0] / EXPECTED['parkes-htru'][1]
norm_rate = norm_sim_rate / norm_real_rate

for i, surv_pop in enumerate(surv_pops):

    name = surv_pop.name.split('_')[-1]
    pprint(name)
    if surv_pop.n_sources() == 0:
        print(surv_pop.source_rate)
        print(f'{name} | no FRBs in population')
        continue

    names.append(name)
    ys.append(y)

    # Dimensions measure plot
    ax = axes[0]
    ax.set_xlabel(r'DM ($\textrm{pc}\ \textrm{cm}^{-3}$)')
    ax.set_ylabel(r'\#')
    ax.set_yscale('log')

    bins, values = hist(surv_pop.frbs.dm, bin_type='lin', norm='frac',
                        n_bins=20)
    values = values.astype(np.float64)
    values *= float(surv_pop.source_rate.f_area)*1e6
    ax.step(bins, values, where='mid', label=name)

    # Fluence plot
    ax = axes[1]
    ax.set_xlabel('S/N')
    ax.set_xscale('log')
    ax.set_ylabel(r'\#(${>}\text{S/N}$)')
    ax.set_yscale('log')

    # Update fluence plot
    bins, values = hist(surv_pop.frbs.snr, bin_type='log', norm='frac',
                        n_bins=25)

    # Cumulative sum
    values = np.cumsum(values[::-1])[::-1]
    values = values.astype(np.float64)
    values *= float(surv_pop.source_rate.f_area)*1e6
    ax.step(bins, values, where='mid', label=name)

    # Plot rates
    ax = axes[2]
    ax.set_xscale('log')
    ax.set_xlabel(r'Rate (day$^{-1}$)')
    rate = surv_pop.source_rate.det/norm_rate
    print(f'rate: {rate}')
    line = ax.errorbar(rate, y,
                       fmt='x',
                       label=rf'{name}')
    ax.grid()
    rates.append(rate)

    y += 1

ax.yaxis.tick_right()
ax.set_yticks(ys)
ax.set_yticklabels(names)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
for i, y in enumerate(ax.get_yticklabels()):
    y.set_color(colors[i])
ax.invert_yaxis()  # labels read top-to-bottom

# Add thin grey horizontal lines
x_lim = ax.get_xlim()
ax.set_xlim(x_lim)
for i, y in enumerate(ys):
    ax.plot((x_lim[0], rates[i]), (y, y), color='k', lw=0.5, zorder=0, ls='--')

for e in list(zip(SURVEYS, rates)):
    pprint(e)

euclidean_lines = True
if euclidean_lines:
    xlims = axes[1].get_xlim()
    ylims = axes[1].get_ylim()
    axes[1].set_xlim(xlims)
    axes[1].set_ylim(ylims)
    xs = np.logspace(np.log10(xlims[0]),
                     np.log10(xlims[1]),
                     100)
    for n in range(-10, 15):
        ys = 10**((np.log10(xs)+n)*-1.5)
        axes[1].plot(xs, ys, 'k:', linewidth=0.25)

# plt.legend()
plt.tight_layout()
plt.savefig(rel_path('./plots/future_surveys.pdf'))
