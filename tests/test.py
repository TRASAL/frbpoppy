import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(precision=2)


def test_other(n, simple=False):

    # Each column is a burst, each row a source
    if simple:
        m = 5
        times = np.random.random((n, m))
        max_time = 1
    else:
        r = 5.7
        k = 0.34
        m = int(round(np.log10(n))*2-1)
        if m < 1:
            m = 1
        print(n, m)
        # Ensure intra-column dependance
        times = 86400*r*np.random.weibull(k, (n, m)).astype(np.float32)
        max_time = 12*60*60
        plt.yscale('log')

    times = np.cumsum(times,axis=1)
    cond = times < max_time
    times[~cond] = np.nan
    bursts = ~np.isnan(times)
    n_bursts = np.count_nonzero(bursts, axis=1)
    unique, counts = np.unique(n_bursts, return_counts=True)
    print(dict(zip(unique, counts)))
    plt.hist(n_bursts, bins=np.linspace(0, m, m+1), density=True)
    plt.show()

test_other(int(1e6))
