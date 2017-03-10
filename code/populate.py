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
             cosmology=True,
             cosmo_pars=[69.6, 0.286, 0.714],
             z_max=8.0,
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
    pop.name = 'Initial'
    pop.n_gen = n_gen
    pop.electron_model = electron_model
    pop.cosmology = cosmology
    pop.H_0 = cosmo_pars[0]
    pop.W_m = cosmo_pars[1]
    pop.W_v = cosmo_pars[2]
    pop.z_max = z_max
    pop.v_max = go.z_to_v(z_max)
    pop.lum_min = lum_dist_pars[0]
    pop.lum_max = lum_dist_pars[1]
    pop.lum_pow = lum_dist_pars[2]
    pop.si_mean = si_pars[0]
    pop.si_sigma = si_pars[1]

    # Create a comoving distance to redshift lookup table
    ds, zs = go.dist_lookup(cosmology=pop.cosmology,
                            H_0=pop.H_0,
                            W_m=pop.W_m,
                            W_v=pop.W_v,
                            z_max=pop.z_max
                            )

    while pop.n_srcs < pop.n_gen:

        # Initialise
        src = Source()

        # Add random coordinates
        src.gb = math.degrees(math.asin(random.random()))
        if random.random() < 0.5:
            src.gb *= -1
        src.gl = random.random() * 360.0

        # Convert coordinates
        # Calculate comoving distance [Gpc]
        src.dist = (pop.v_max * random.random() * (3/(4*math.pi)))**(1/3)
        src.z = go.interpolate_z(src.dist, ds, zs, H_0=pop.H_0)
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
        src.w_int = (1+src.z)*random.uniform(0.1,10)

        # Add bolometric luminosity [W]
        src.lum_bol = 1e64  # 8.0e37
        #src.lum_bol = ds.powerlaw(pop.lum_min, pop.lum_max, pop.lum_pow)

        # Add spectral index
        src.si = -1.4
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
