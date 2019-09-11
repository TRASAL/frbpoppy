from collections import defaultdict
import numpy as np

# Everything in fractions of a day
times = np.array([0.1])  # .astype(np.float32)
frac_vis = 0.5

def cal_p_visible(times, frac):
    p_det = defaultdict(float)

    lower = times - frac
    lower[(lower < 0)] = 0.
    upper = times + frac
    upper[(upper > 1.)] = 1.

    # Determine times at which something changes
    t_change = np.hstack((lower, upper))
    # What changes in number of bursts seen
    change = np.hstack((np.ones_like(times), -1*np.ones_like(times)))

    # Sort times
    s = t_change.argsort()
    t_change = t_change[s]
    change = change[s]

    # Fraction of time in a certain state
    frac_change = np.diff(t_change)
    state_change = np.cumsum(change)[:-1]

    # Add together all fractions of beams
    mapping = {}
    for i in np.unique(state_change):
        if i != 0:
            p_det[int(i)] += sum(frac_change[np.where(state_change == i)[0]])

    p_det[0] = 1. - sum(p_det.values())

    return p_det


if __name__ == '__main__':
    import unittest

    p = cal_p_visible

    class TestVisibilities(unittest.TestCase):

        def test_single(self):
            times = np.array([0.5])
            ps = p(times, 0.5)
            self.assertAlmostEqual(ps[0], 0)
            self.assertAlmostEqual(ps[1], 1)

        def test_boundary(self):
            times = np.array([0.1])
            ps = p(times, 0.5)
            self.assertAlmostEqual(ps[0], 0.4)
            self.assertAlmostEqual(ps[1], 0.6)

        def test_two(self):
            times = np.array([0.2, 0.8])
            ps = p(times, 0.1)
            self.assertAlmostEqual(ps[0], 0.6)
            self.assertAlmostEqual(ps[1], 0.4)

        def test_two_close(self):
            times = np.array([0.4, 0.6])
            ps = p(times, 0.3)
            self.assertAlmostEqual(ps[0], 0.2)
            self.assertAlmostEqual(ps[1], 0.4)
            self.assertAlmostEqual(ps[2], 0.4)

        def test_three_close(self):
            times = np.array([0.4, 0.5, 0.6])
            ps = p(times, 0.3)
            self.assertAlmostEqual(ps[0], 0.2)
            self.assertAlmostEqual(ps[1], 0.2)
            self.assertAlmostEqual(ps[2], 0.2)
            self.assertAlmostEqual(ps[3], 0.4)

        def test_three_close_boundary(self):
            times = np.array([0.7, 0.8, 0.9])
            ps = p(times, 0.3)
            self.assertAlmostEqual(ps[0], 0.4)
            self.assertAlmostEqual(ps[1], 0.1)
            self.assertAlmostEqual(ps[2], 0.1)
            self.assertAlmostEqual(ps[3], 0.4)


    unittest.main(verbosity=2)
