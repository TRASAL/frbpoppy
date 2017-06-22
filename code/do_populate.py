"""Code for generating a population of FRBs."""

import math
import random

import galacticops as go
import distributions as dis
import precalc as pc
from population import Population
from source import Source


def generate(n_gen,
             days=1,
             cosmology=True,
             cosmo_pars=[69.6, 0.286, 0.714],
             electron_model='ne2001',
             emission_pars=[10e6, 10e9],
             lum_dist_pars=[1e40, 1e50, 1],
             name=None,
             pulse=[0.1, 5],
             repeat=0.0,
             si_pars=[-1.4, 0.],
             z_max=2.5,
             test=False):

    """Generate a population of FRB sources

    Args:
        n_gen (int): Number of FRB sources/sky/time to generate
        days (float): Number of days over which FRBs are generated.
            Defaults to 1
        cosmology (bool, optional): Whether to use cosmology in all
            calculations. Defaults to True.
        cosmo_pars (list, optional): Three values, the first being the Hubble
            constant, the second the density parameter Ω_m and the third
            being the cosmological constant Ω_Lambda (referred to as W_m
            in the rest of the code). These parameters default to those
            determined with Planck [69.6, 0.286, 0.714]
        electron_model (str, optional): Model for the free electron density in
            the Milky Way. Defaults to 'ne2001'
        emission_pars (list, optional): The minimum and maximum frequency [Hz]
            between which FRB sources should emit the given bolometric
            luminosity. Defaults to [10e6,10e9]
        lum_dist_pars (list, optional): Bolometric luminosity distribution
            parameters being: the minimum luminosity [erg/s], the maximum
            luminosity [erg/s] and the powerlaw index of the distribution.
            Defaults to [1e50, 1e90, 1]
        name (str, optional): Name to be given to the population. Defaults to
            'initial'
        pulse (str, optional): Values between which the intrinsic pulse width
            should be [ms]. Defaults to [0.5, 5]
        repeat (float, optional): Fraction of sources which repeat
        si_pars (list, optional): Spectral index parameters, being the mean
            index and the standard deviation thereof. Defaults to [-1.4, 0.]
        z_max (float, optional): The maximum redshift out to which to
            distribute FRBs
        test (float, optional): Flag to help testing

    Returns:
        pop (Population): A population of generated sources

    """

    # Check input
    if type(n_gen) != int:
        m = 'Please ensure the number of generated sources is an integer'
        raise ValueError(m)

    if n_gen < 1:
        m = 'Please ensure the number of generated sources is > 0'
        raise ValueError(m)

    # Check input
    if type(days) != int:
        m = 'Please ensure the number of days is an integer'
        raise ValueError(m)

    if type(cosmology) != bool:
        raise ValueError('Please ensure cosmology is a boolean')

    if not all(isinstance(par, (float, int)) for par in cosmo_pars):
        m = 'Please ensure all cosmology parameters are floats or integeters'
        raise ValueError(m)

    if len(cosmo_pars) != 3:
        m = 'Please ensure there are three cosmology parameters'
        raise ValueError(m)

    if electron_model not in ['ne2001']:
        m = 'Unsupported electron model: {}'.format(electron_model)
        raise ValueError(m)

    if not all(isinstance(par, (float, int)) for par in emission_pars):
        m = 'Please ensure all emission parameters are floats or integeters'
        raise ValueError(m)

    if len(emission_pars) != 2:
        m = 'Please ensure there are two emission parameters'
        raise ValueError(m)

    if not all(isinstance(par, (float, int)) for par in lum_dist_pars):
        m = 'Please ensure all luminosity distribution parameters are '
        m += 'floats or integeters'
        raise ValueError(m)

    if len(lum_dist_pars) != 3:
        m = 'Please ensure there are three luminosity distribution parameters'
        raise ValueError(m)

    if not all(isinstance(par, (float, int)) for par in pulse):
        m = 'Please ensure all pulse parameters are '
        m += 'floats or integeters'
        raise ValueError(m)

    if len(pulse) != 2:
        m = 'Please ensure there are two pulse parameters'
        raise ValueError(m)

    if not name:
        name = 'Initial'

    if type(name) != str:
        m = 'Please provide a string as input for the name of the population'
        raise ValueError(m)

    if type(repeat) != float:
        m = 'Please ensure fraction of sources which repeat is a float'
        raise ValueError(m)

    if repeat > 1:
        m = 'The repeat fraction can not be more than 1.0'
        raise ValueError(m)

    if not all(isinstance(par, (float, int)) for par in si_pars):
        m = 'Please ensure all spectral index parameters are '
        m += 'floats or integers'
        raise ValueError(m)

    if len(si_pars) != 2:
        m = 'Please ensure there are two spectral index parameters'
        raise ValueError(m)

    if not isinstance(z_max, (float, int)):
        m = 'Please ensure the maximum redshift is given as a float or integer'
        raise ValueError(m)

    # Set up population
    pop = Population()
    pop.cosmology = cosmology
    pop.electron_model = electron_model
    pop.f_max = emission_pars[1]
    pop.f_min = emission_pars[0]
    pop.H_0 = cosmo_pars[0]
    pop.lum_max = lum_dist_pars[1]
    pop.lum_min = lum_dist_pars[0]
    pop.lum_pow = lum_dist_pars[2]
    pop.name = name
    pop.n_gen = n_gen
    pop.repeat = repeat
    pop.si_mean = si_pars[0]
    pop.si_sigma = si_pars[1]
    pop.time = days * 86400  # Convert to seconds
    pop.v_max = go.z_to_v(z_max)
    pop.w_max = pulse[1]
    pop.W_m = cosmo_pars[1]
    pop.w_min = pulse[0]
    pop.W_v = cosmo_pars[2]
    pop.z_max = z_max

    while pop.n_srcs < pop.n_gen:

        # Initialise
        src = Source()

        # Add random directional coordinates
        src.gl = random.random() * 360.0 - 180
        src.gb = math.degrees(math.asin(random.random()))
        if random.random() < 0.5:
            src.gb *= -1

        # Convert
        src.ra, src.dec = go.lb_to_radec(src.gl, src.gb)

        # Calculate comoving distance [Gpc]
        src.dist = (pop.v_max * random.random() * (3/(4*math.pi)))**(1/3)

        # Calculate redshift
        src.z = pc.dist_table(src.dist, H_0=pop.H_0, test=test)

        # Convert into galactic coordinates
        src.gx, src.gy, src.gz = go.lb_to_xyz(src.gl, src.gb, src.dist)

        # Dispersion measure of the Milky Way
        src.dm_mw = pc.ne2001_table(src.gl, src.gb, test=test)

        # Dispersion measure of the intergalactic medium
        src.dm_igm = go.ioka_dm_igm(src.z)

        # Dispersion measure of the host
        src.dm_host = 100.  # Thornton et al. (2013)

        # Total dispersion measure
        src.dm = src.dm_mw + src.dm_igm + src.dm_host

        # Add initial frb
        src.create_frb(pop)

        # If repeating add another FRB
        if random.random() < pop.repeat:

            ts = dis.oppermann_pen()

            for t in ts:
                src.create_frb(pop, time=t)

        # Add source to population
        pop.add(src)

    # Save population
    pop.pickle_pop()

    return pop
