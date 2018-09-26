import numpy as np
from bokeh.plotting import figure, show
from frbpoppy import Frbcat, Survey

frbcat = Frbcat().df
frbcat = frbcat[(frbcat.survey == 'HTRU')]
pop = Frbcat().to_pop(df=frbcat)

w_effs = pop.get('w_eff')
dms = pop.get('dm')

survey = Survey('HTRU')
t_dms = [8.297616e6*survey.bw_chan*dm*(survey.central_freq)**-3 for dm in dms]

ratios = [e[0]/e[1] for e in zip(w_effs, t_dms)]

# Bin up
hist, edges = np.histogram(ratios, bins=20)
cum_hist = [sum(hist[i:]) for i in range(len(hist))]

# Plot
p = figure(title='Ratio tau / tau_dm', x_axis_type='linear', y_axis_type='log')
p.yaxis.axis_label = 'Fraction'

# Find optimum lower limit
m = min(e for e in cum_hist if e > 0)
bottom = 10**(np.log10(m) - 1)

p.quad(top=hist,
       bottom=bottom,
       left=edges[:-1],
       right=edges[1:],
       alpha=0.5,
       legend='HTRU')

show(p)
