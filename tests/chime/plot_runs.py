"""Plot information on various CHIME runs"""
import pandas as pd
import matplotlib.pyplot as plt

from tests.convenience import rel_path

df = pd.read_csv(rel_path('chime/runs.csv'), index_col=0)

exp_one = 200
exp_rep = 2.5

fig, axes = plt.subplots(1, 2, sharex=True, sharey=True)
ax1, ax2 = axes

max_clb = 200
c = ax1.scatter(df.srcs, df.rate, c=df.one_offs-exp_one, s=500, cmap='RdBu',
                vmin=-max_clb, vmax=max_clb)
fig.colorbar(c, ax=ax1)
ax1.set_xscale('log')
ax1.set_xlabel('# sources')
ax1.set_title('One-offs')
ax1.set_ylabel('rate (/day)')

max_clb = 50
d = ax2.scatter(df.srcs, df.rate, c=df.repeaters-exp_rep, s=500, cmap='RdBu',
                vmin=-30, vmax=30)
ax2.set_xscale('log')
ax2.set_xlabel('# sources')
ax2.set_title('Repeaters')

for i, txt in enumerate(df.index):
    ax1.annotate(txt, (df.srcs.iat[i], df.rate.iat[i]))
    ax2.annotate(txt, (df.srcs.iat[i], df.rate.iat[i]))

clb = fig.colorbar(d, ax=ax2, label='Differene w.r.t. sought value')
plt.show()

# Plot change in rates between multiple generations
db = df[(df.srcs == 34000) & (df.rate == 9.0)]
fig, axes = plt.subplots(1, 2, sharey=True)
ax1, ax2 = axes
ax1.hist(db.one_offs)
ax1.set_xlabel('One-offs')
ax2.hist(db.repeaters)
ax2.set_xlabel('Repeaters')
plt.show()
