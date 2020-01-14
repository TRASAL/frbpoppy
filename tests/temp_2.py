import numpy as np
import bisect
import time

# Initialize an example of an array in which to search
r, c = int(1e2), int(1e6)
a = np.arange(r*c).reshape(r, c)

# Set up search limits
min_v = float(1e24)
max_v = float(1e24)


# Find indices of occurances
def use_np_where():
    return np.where(((a >= min_v) & (a <= max_v)))


def use_searchsorted():
    i1 = np.searchsorted(a.ravel(), min_v, 'left')
    i2 = np.searchsorted(a.ravel(), max_v, 'right')
    return np.unravel_index(np.arange(i1, i2), a.shape)


def use_bisect():
    import IPython; IPython.embed()
    i1 = bisect.bisect_left(a.base, min_v)
    i2 = bisect.bisect_right(a.base, max_v)
    return np.unravel_index(np.arange(i1, i2), a.shape)


def timeit(f):
    print(f.__name__)
    start = time.time()
    i = f()
    end = time.time()
    print(end - start)
    return i


idx = timeit(use_np_where)
idx1 = timeit(use_searchsorted)
idx2 = timeit(use_bisect)
print((idx[0] == idx2[0]).all() and (idx[1] == idx2[1]).all())
