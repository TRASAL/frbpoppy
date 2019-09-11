import numpy as np
from frbpoppy.galacticops import separation


def test_angsep():
    """
    Tests that the angular separation behaves correctly.
    """
    # lon1, lat1, lon2, lat2 in degrees
    coords = [(1, 0, 0, 0),
              (0, 1, 0, 0),
              (0, 0, 1, 0),
              (0, 0, 0, 1),
              (0, 0, 10, 0),
              (0, 0, 90, 0),
              (0, 0, 180, 0),
              (0, 45, 0, -45),
              (0, 60, 0, -30),
              (-135, -15, 45, 15),
              (100, -89, -80, 89),
              (0, 0, 0, 0),
              (0, 0, 1. / 60., 1. / 60.)]
    correct_seps = [1, 1, 1, 1, 10, 90, 180, 90, 90, 180, 180, 0,
                    0.023570225877234643]
    correctness_margin = 2e-10

    for (lon1, lat1, lon2, lat2), corrsep in zip(coords, correct_seps):
        angsep = separation(lon1, lat1, lon2, lat2)
        print(lon1, lat1, lon2, lat2, f'Should be {corrsep}, is {angsep}')
        assert np.fabs(angsep - corrsep) < correctness_margin


test_angsep()
