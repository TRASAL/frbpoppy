"""Code for generating a population of FRBs"""

import math
import random
import sys
from argparse import ArgumentParser

from log import Log
from population import Population
from source import Source

assert sys.version_info >= (3,0), 'Please run with Python3'


def generate(n_gen,
             log_loc=None,
             no_log=False,
             quiet=False,
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

    # Log parameters
    pop.logloc = log_loc
    pop.no_log = no_log
    pop.quiet = quiet
    pop.verbose = verbose

    while pop.n_det < pop.n_gen:

        # Initialise
        src = Source()
        if pop.n_det < 5:
            # Add random coordinates
            src.gb = math.degrees(math.asin(random.random()))
            if random.random() < 0.5:
                src.gb *= -1
            src.gl = random.random() * 360.0

            # Convert coordinates
            src.lb_to_xyz(1.0) #1 kpc

        # Add to population
        pop.population.append(src)
        pop.n_det += 1

    pop.write_csv('logs/test.csv')


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
