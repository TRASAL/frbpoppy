"""Code for generating a population of FRBs"""

import math
import random
from argparse import ArgumentParser

import distributions as ds
import galacticops as go
from log import pprint
from population import Population
from source import Source


def generate(n_gen,
             electron_model='ne2001',
             z_max=2.5,
             lum_dist_pars=[1, 2, 1],
             si_pars=[1,1],
             scindex=-3.86):

    """
    Generate a population of FRB sources

    Args:
        ngen (int): Number of FRB sources to generate
    """
    # Check input

    pop = Population()
    pop.n_gen = n_gen
    pop.lum_dist_pars = lum_dist_pars
    pop.si_pars = si_pars

    pop.electron_model = electron_model

    while pop.n_srcs < pop.n_gen:

        # Initialise
        src = Source()

        # Add random coordinates
        src.gb = math.degrees(math.asin(random.random()))
        if random.random() < 0.5:
            src.gb *= -1
        src.gl = random.random() * 360.0

        # Convert coordinates
        src.z = z_max*random.random()
        src.dist = go.z_to_d(src.z)  # [kpc]
        src.gx, src.gy, src.gz = go.lb_to_xyz(src.gl, src.gb, src.dist)

        # Calculate dispersion measure
        # Milky Way
        src.dm_mw = go.ne2001_dist_to_dm(src.dist, src.gl, src.gb)
        # Intergalactic medium
        src.dm_igm = go.ioka_dm_igm(src.z)
        # Host
        src.dm_host = 100.  # Thornton et al. (2013)
        # Total
        src.dm = src.dm_mw + src.dm_igm + src.dm_host

        # Give an intrinsic pulse width [ms]
        src.w_int = random.uniform(0.1,10)

        # Add luminosity at 1400 MHz
        src.lum_1400 = 10.0#1.14
        #src.lum_1400 = ds.powerlaw(pop.lum_min, pop.lum_max, pop.lum_pow)

        # Add spectral index
        src.si = 1.0
        #src.si = random.gauss(pop.si_mean, pop.si_sigma)

        # Add to population
        pop.sources.append(src)
        pop.n_srcs += 1

    return pop


if __name__ == '__main__':

    parser = ArgumentParser(description='Generate a population of FRB sources')

    # Scientific parameters
    parser.add_argument('-n',
                        '--n_gen',
                        type=int,
                        required=False,
                        help='number of FRB sources to generate/detect')

    # Logging options
    parser.add_argument('-nl',
                        '--no_log',
                        action='store_true',
                        help="don't save log to file")

    parser.add_argument('-ll',
                        '--log_loc',
                        default=None,
                        help='location of log file')

    # Verbosity
    parser.add_argument('-v',
                        '--verbose',
                        action='store_true',
                        help='get more output in the console')

    parser.add_argument('-q',
                        '--quiet',
                        action='store_true',
                        help='turn off output in the console')

    args = parser.parse_args()

    generate(n_gen=args.n_gen,
             no_log=args.no_log,
             verbose=args.verbose,
             quiet=args.quiet,
             log_loc=args.log_loc)
