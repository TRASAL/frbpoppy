"""Monte Carlo over various frbpoppy parameters."""

from collections import namedtuple
import numpy as np
import pandas as pd

from do_populate import generate
from do_survey import observe

# If set is None, this parameter refers to a range
Par = namedtuple('par', 'min max step set')

dm_host = Par(0, 200, 10, 100)
dm_igm_slope = Par(1000, 1400, 20, 1200)
freq_min = Par(10e5, 10e10, 0.5e1, 10e6)
freq_max = Par(10e5, 10e10, 0.5e1, 10e9)
lum_bol_slope = Par(0.5, 1.5, 0.1, 1.)
lum_bol_min = Par(1e30, 1e60, 1e1, 1e40)
lum_bol_max = Par(1e30, 1e60, 1e1, 1e50)
n_day = Par(2000, 14000, 2000, 10000)
repeat = Par(0, 0.1, 0.01, 0.05)
si_mean = Par(-2.0, -1, 0.1, -1.4)
si_sigma = Par(0.0, 0.5, 0.1, 0.0)
w_int_min = Par(0.1, 5, 0.1, 1.)
w_int_max = Par(0.1, 5, 0.1, 5.)

pars = {'dm_host': dm_host,
        'dm_igm_slope': dm_igm_slope,
        'freq_min': freq_min,
        'freq_max': freq_max,
        'lum_bol_slope': lum_bol_slope,
        'lum_bol_min': lum_bol_min,
        'lum_bol_max': lum_bol_max,
        'n_day': n_day,
        'repeat': repeat,
        'si_mean': si_mean,
        'si_sigma': si_sigma,
        'w_int_min': w_int_min,
        'w_int_max': w_int_max}


def par_range(p, mi=None, ma=None, st=None, se=None):
    """Quick range generator."""
    if not mi:
        mi = p.min
    if not ma:
        ma = p.max
    if not st:
        st = p.step
    if not se:
        se = p.set
    return np.arange(mi, ma+st, st)

# Set all values
inputs = {k:[] for k in pars}
print(inputs)

# Loop through combinations of min and max pulse widths
for mi in par_range(w_int_min):
    for ma in par_range(w_int_max, mi=mi):

        for p in pars:
            if p == 'w_int_min':
                inputs[p].append(mi)
            elif p == 'w_int_max':
                inputs[p].append(ma)
            else:
                inputs[p].append(pars[p].set)

df = pd.DataFrame(inputs)
print(df)
