import os
import unittest

import distributions as ds
import galacticops as go
from population import Population as pop
from do_populate import generate
from do_survey import observe
from do_plot import plot


class TestFullRun(unittest.TestCase):

    def test_run(self):
        # Test a full run of frbpoppy
        # Should be the equivalent to a bucket load of rainbows, unicorns and
        # puppies. All other tests just check the weird cases work.
        population = generate(10)
        survey_population = observe(population, 'WHOLESKY')
        plot(population, survey_population, show=False)


class TestDistributions(unittest.TestCase):

    def test_powerlaw(self):

        # Check normal usage
        r1 = ds.powerlaw(1, 2, 2)
        self.assertTrue(1 <= r1 <= 2)

        # Check inverted input
        r2 = ds.powerlaw(2, 1, 2)
        self.assertTrue(1 <= r1 <= 2)

        # Check what happens when zero is raised to a negative power
        with self.assertRaises(ValueError):
            ds.powerlaw(0, 5, -2)


class TestPopulate(unittest.TestCase):

    def testgenerate(self):

        # Check generating a population works
        r1 = generate(2, name='temp')
        self.assertIsInstance(r1, pop)

        # Check printing a population works
        self.assertTrue(r1.__str__())

        # Check saving works
        r1.save()

        # Check different input values
        with self.assertRaises(ValueError):
            generate(0)
        with self.assertRaises(ValueError):
            generate(1.5)
        with self.assertRaises(ValueError):
            generate(1, 'string')
        with self.assertRaises(ValueError):
            generate(1, cosmology=1)
        with self.assertRaises(ValueError):
            generate(1, cosmo_pars=[1, 2])
        with self.assertRaises(ValueError):
            generate(1, cosmo_pars=[True, '1'])
        with self.assertRaises(ValueError):
            generate(1, emission_pars=[2])
        with self.assertRaises(ValueError):
            generate(1, emission_pars=[2, 'test'])
        with self.assertRaises(ValueError):
            generate(1, lum_dist_pars=[True, '1'])
        with self.assertRaises(ValueError):
            generate(1, lum_dist_pars=[1, 2])
        with self.assertRaises(ValueError):
            generate(1, name=1.5)
        with self.assertRaises(ValueError):
            generate(1, si_pars=[True])
        with self.assertRaises(ValueError):
            generate(1, si_pars=[1.2, True, 'test'])
        with self.assertRaises(ValueError):
            generate(1, z_max=[0, 0])
        with self.assertRaises(ValueError):
            generate(1, electron_model='unsupported')

    def tearDown(self):
        """Clean up temp files"""
        loc = '../data/results/population_temp.csv'
        out = os.path.join(os.path.dirname(__file__), loc)
        os.remove(out)


class TestSurvey(unittest.TestCase):

    def setUp(self):

        self.pop = generate(10)
        self.bright_pop = generate(10, lum_dist_pars=[1e80, 1e90, 1])

    def testobserve(self):

        # Check different input values
        with self.assertRaises(IOError):
            observe(self.pop, 'ABCDEFG')

        # Check printing a survey works
        s = observe(self.pop, 'PMSURV', return_pop=False)
        self.assertTrue(s.__str__())

        # Check conducting a survey works
        s1 = observe(self.pop, 'PMSURV')
        self.assertIsInstance(s1, pop)

        # Check bright FRBs are always found
        s2 = observe(self.bright_pop, 'WHOLESKY')
        self.assertEqual(len(s2.sources), len(self.pop.sources))

        # Not the best tests for scat and scint, but will have to do
        # Check scattering works
        s3 = observe(self.bright_pop, 'WHOLESKY', scat=True)
        self.assertIsInstance(s3, pop)

        # Check scintillation works
        s4 = observe(self.bright_pop, 'WHOLESKY', scint=True)
        self.assertIsInstance(s4, pop)

        # Test printing of sources
        self.assertTrue(s3.sources[0].__str__())


class TestGalacticops(unittest.TestCase):

    def test_radec_to_lb(self):
        """Test values from Robert Martin Ayers's calculator"""

        # Actual FRB coordinates
        gl, gb = go.radec_to_lb('19:06:53', '-40:37:14')
        self.assertTrue(abs(gl + (360 - 356.6411)) <= 1.5)
        self.assertTrue(abs(gb + 20.0207) <= 1.5)

        # Some extreme cases
        gl, gb = go.radec_to_lb('0:0:0', '0:0:0')
        self.assertTrue(abs(gl - 96.3373) <= 1.5)
        self.assertTrue(abs(gb + 60.1886) <= 1.5)

        gl, gb = go.radec_to_lb('23:0:0', '-90:0:0')
        self.assertTrue(abs(gl + (360 - 302.9319)) <= 1.5)
        self.assertTrue(abs(gb + 27.1283) <= 1.5)


if __name__ == '__main__':
    unittest.main()
