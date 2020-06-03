# Plot CHIME runs
import pandas as pd
import mpl_toolkits.mplot3d
import matplotlib.pyplot as plt

from tests.convenience import rel_path

df = pd.read_csv(rel_path('chime_runs.csv'), index_col=0)

exp_one = 200
exp_rep = 2.5

# fig = plt.figure()
# ax1 = fig.add_subplot(121, projection='3d')
#
# cmap = plt.cm.get_cmap('viridis')
#
# ax1.scatter(df.srcs,
#             df.rate,
#             df.one_offs,
#             label='One Offs',
#             marker='x',
#             s=50,
#             cmap=cmap,
#             c=range(len(df)))
#
# ax1.set_xlabel('# sources')
# ax1.set_ylabel('rate (/day)')
# ax1.set_zlabel('# one-off detections')
#
# ax2 = fig.add_subplot(122, projection='3d')
#
# c = ax2.scatter(df.srcs,
#                 df.rate,
#                 df.repeaters,
#                 label='Repeaters',
#                 marker='o',
#                 s=50,
#                 cmap=cmap,
#                 c=range(len(df)))
#
#
# ax2.set_xlabel('# sources')
# ax2.set_ylabel('rate (/day)')
# ax2.set_zlabel('# repeater detections')
#
# plt.colorbar(c)
# plt.show()

# fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True)
#
# ax1.scatter(df.index, (df.repeaters/exp_rep), label='repeaters')
# ax1.scatter(df.index, (df.one_offs/exp_one), label='one-offs')
# ax1.set_ylabel('Fraction of expectation value')
#
# ax2.scatter(df.index, (df.repeaters-exp_rep), label='repeaters')
# ax2.scatter(df.index, (df.one_offs-exp_one), label='one-offs')
# ax2.set_ylabel('Difference with expectation value')
#
# ax3.scatter(df.index, df.one_offs, label='one-offs')
# ax3.axhline(200, c='r', ls='-')
# ax3.set_ylabel('Value')
#
# ax4.scatter(df.index, [0]*len(df.index), label='one-offs', alpha=0)
# ax4.scatter(df.index, df.repeaters, label='repeaters')
# ax4.axhline(2.5, c='r', ls='-')
# ax4.set_ylabel('Value')
#
# ax4.set_xlabel('Run number')
# ax1.legend()
# plt.show()

fig, axes = plt.subplots(1, 2, sharey=True)
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
# clb.ax.set_title('Difference w.r.t. sought value')
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
