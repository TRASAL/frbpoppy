"""Code for generating a population of FRBs"""

import math
import random

from argparse import ArgumentParser

import distributions as ds
import galacticops as go
from population import Population
from source import Source


def generate(n_gen,
             electron_model='ne2001',
             log_loc=None,
             lum_dist_pars=[0, 0, 0],
             no_log=False,
             quiet=False,
             scindex=-3.86,
             verbose=False):

    """
    Generate a population of FRB sources

    Args:
        ngen (int): Number of FRB sources to generate
        no_log (boolean): Don't save a log of this run
        log_loc (str): Location of log filename
        verbose (boolean): Up the verbosity
        quiet (boolean): Turn off piping output to console
    """

    pop = Population()

    pop.n_gen = n_gen
    pop.lum_dist_pars = lum_dist_pars

    # Log parameters
    pop.logloc = log_loc
    pop.no_log = no_log
    pop.quiet = quiet
    pop.verbose = verbose

    while pop.n_det < pop.n_gen:

        # Initialise
        src = Source()

        # Add random coordinates
        src.gb = math.degrees(math.asin(random.random()))
        if random.random() < 0.5:
            src.gb *= -1
        src.gl = random.random() * 360.0

        # Convert coordinates
        d = 5*random.random()  # [kpc]
        src.gx, src.gy, src.gz = go.lb_to_xyz(src.gl, src.gb, d)

        # Calculate distance to source
        src.dist = go.calc_d_sun(src.gx, src.gy, src.gz)

        # Calculate dispersion measure
        src.dm = go.ne2001_dist_to_dm(src.dist, src.gl, src.gb)

        # Calculate intrinsic pulse width [ms]
        src.width = 3.0

        # Add luminosity at 1400 MHz
        src.lum_1400 = 1.14
        #src.lum_1400 = ds.powerlaw(pop.lum_min, pop.lum_max, pop.lum_pow)

        # Add to population
        pop.sources.append(src)
        pop.n_det += 1

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
