"""Plot various parameter distribution for different beam types."""
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from frbpoppy.do_populate import generate
from frbpoppy.do_survey import observe
from frbpoppy.population import unpickle

MAKE = False  # Construct a population to survey
PARMS = ['fluence']  # , 'fluence', 's_peak', 'w_eff']
PATTERNS = ['perfect', 'tophat', 'gaussian', 'airy']
SURVEY = 'HTRU'
NBINS = 50

if MAKE:
    n_per_day = 5000
    days = 28
    pop_std = generate(n_per_day*days, days=days, name='standard')
else:
    pop_std = unpickle('standard')


for p in PARMS:

    fig = plt.figure()
    ax = fig.add_subplot(111)

    if p == 'fluence':
        title = r'Fluence (Jy ms)'
        plt.xscale('log')
        # plt.yscale('log')
        logbins = True
    elif p == 's_peak':
        title = r'S_{peak}'
        plt.xscale('log')
        plt.yscale('log')
        logbins = True
    elif p == 'dm':
        title = 'Dispersion Measure'
    elif p == 'w_eff':
        title = 'Pulse Width'
    else:
        logbins = False

    for pattern in PATTERNS:
        surv_pop = observe(pop_std, SURVEY, gain_pattern=pattern)

        ps = []

        for src in surv_pop.sources:
            try:
                ps.append(getattr(src, p))
            except AttributeError:
                for frb in src.frbs:
                    ps.append(getattr(frb, p))

        if logbins:
            # Calculate the min and max powers:
            start_power = np.floor(np.log10(min(ps)))
            end_power = np.ceil(np.log10(max(ps)))
            num_bins = (end_power - start_power) + 1
            # Generate a range of delimiters in log space
            bins = np.logspace(start_power, end_power, NBINS, base=10)
        else:
            bins = NBINS

        n, bins, patches = ax.hist(ps, bins=bins, density=False,
                                   label=pattern, histtype='step')

    # Create new legend handles but use the colors from the existing ones
    handles, labels = ax.get_legend_handles_labels()
    new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]

    plt.legend(handles=new_handles, labels=labels)

    plt.xlabel(title)
    plt.ylabel(r'\# FRB Detections')
    plt.tight_layout()
    plt.savefig(f'plots/{p}_per_beampattern.pdf')
