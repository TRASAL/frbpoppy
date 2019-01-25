import numpy as np
import matplotlib.pyplot as plt

from frbpoppy import distributions

def plot(power):
    plt.clf()
    pl = distributions.powerlaw(1e35, 1e40, power, 100000)
    minx = min(pl)
    maxx = max(pl)
    bins = 10 ** np.linspace(np.log10(minx), np.log10(maxx), 50)
    plt.hist(pl, bins=bins)
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig('plots/powerlaw.pdf')


plot(-0.05)
