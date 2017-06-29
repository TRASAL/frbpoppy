"""Run a Monte Carlo simulation over a wide range of parameters."""

from collections import defaultdict
from scipy.stats import ks_2samp
import pandas as pd

from do_populate import generate
from do_survey import observe
from frbcat import get_frbcat
from monte_carlo import MonteCarlo

surveys = {'WHOLESKY': 'WHOLESKY',
           'APERTIF': 'APERTIF',
           'PMSURV': 'parkes',
           'HTRU': 'parkes',
           'ASKAP-INCOH': 'ASKAP',
           'ASKAP-FLY': 'ASKAP',
           'GBT': 'GBT',
           'PALFA': 'arecibo',
           'ARECIBO-SPF': 'arecibo',
           'ALFABURST': 'arecibo',
           'UTMOST-1D': 'UTMOST'}

# Get parameters over which to loop
mc = MonteCarlo()
mc.w_int()
mc.save()
df = mc.read().head(5)

# Get the actual observations with which to compare
cat = get_frbcat()

# Set up dictionary for results
d = defaultdict(list)

# Iterate over each set of parameters
for i, r in df.iterrows():

    # Save each set of parameters
    sett = r.to_dict()

    # Create population
    pop = generate(r.n_day,
                   days=1,
                   emission_pars=[r.freq_min, r.freq_max],
                   lum_dist_pars=[r.lum_bol_min,
                                  r.lum_bol_max,
                                  r.lum_bol_slope],
                   pulse=[r.w_int_min, r.w_int_max],
                   repeat=r.rep,
                   si_pars=[r.si_mean, r.si_sigma])

    for s in surveys:
        # Survey population
        sur_pop = observe(pop, s, output=False).to_df()

        # Find matching properties
        cols = [c for c in cat if c in sur_pop]

        # Collect results
        ks = {}

        # KS-test each parameter
        for c in cols:
            obs_cat = cat[(cat.Telescope == surveys[s])][c]
            obs_cat = pd.to_numeric(obs_cat)
            obs_pop = pd.to_numeric(sur_pop[c], errors='coerce')

            ks['ks_' + c] = ks_2samp(obs_pop, obs_cat)[1]

        # Add as results
        for p in sett:
            d[p].append(sett[p])
        for k in ks:
            d[k].append(ks[k])
        d['telescope'].append(surveys[s])
        d['survey'].append(s)

mc.save(df=pd.DataFrame(d), filename='ks.db')
