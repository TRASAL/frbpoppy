"""Test survey strategies."""
import numpy as np
import unittest
from frbpoppy import FRBs, Survey


class TestStrategies(unittest.TestCase):
    """Test survey strategies."""

    def setUp(self):
        """Set up frbs in three pointings."""
        self.frbs = FRBs()

        # Three pointings
        self.frbs.ra = np.array([100, 100, 100])
        self.frbs.dec = np.array([-10, 0, 10])

        self.survey = Survey('chime', n_days=1)
        self.survey.t_obs = int(86400/3)
        self.survey.pointings = list(zip(self.frbs.ra, self.frbs.dec))

    def test_regular(self):
        """Test a regular survey strategy."""
        self.frbs.time = np.array([[.3, .6, 1],  # pointing A
                                   [.3, .6, 1],  # pointing B
                                   [.3, .6, 1]])  # pointing C
        self.frbs.time *= 86400

        answer = np.array([[True, False, False],
                          [False, True, False],
                          [False, False, True]])

        mask = self.survey.in_observation(self.frbs, strategy='regular',
                                          t_stick=self.survey.t_obs)

        np.testing.assert_array_equal(mask, answer)

    def test_followup(self):
        """Test a follow-up survey strategy."""
        self.frbs.time = np.array([[.3, .6, 1],  # pointing A
                                   [.3, .6, 1],  # pointing B
                                   [.3, .6, 1]])  # pointing C
        self.frbs.time *= 86400

        answer = np.array([[True, True, True],
                          [False, False, False],
                          [False, False, False]])

        mask = self.survey.in_observation(self.frbs, strategy='follow-up',
                                          t_stick=self.survey.t_obs)

        np.testing.assert_array_equal(mask, answer)

    def test_followup_with_gap(self):
        """Test a follow-up survey strategy."""
        np.warnings.filterwarnings('ignore')

        self.frbs.time = np.array([[.3, np.nan, 1],  # pointing A
                                   [.3, .6, 1],  # pointing B
                                   [.3, .6, 1]])  # pointing C
        self.frbs.time *= 86400

        answer = np.array([[True, False, False],
                          [False, False, True],
                          [False, False, False]])

        mask = self.survey.in_observation(self.frbs, strategy='follow-up',
                                          t_stick=self.survey.t_obs)

        np.testing.assert_array_equal(mask, answer)


if __name__ == '__main__':
    unittest.main()
