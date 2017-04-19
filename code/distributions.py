import random


def powerlaw(low, high, power):
    """
    Return random variable distributed according to power law

    Args:
        low (float): Lower limit of distribution
        high (float): Higher limit of distribution
        power (float): Power of power law distribution
    Returns:
        x (float): Random variable picked from power law distribution
    """
    if low > high:
        low, high = high, low

    p1 = power + 1
    if (low == 0. or high == 0.) and p1 < 0:
        raise ValueError('Power laws are not defined at 0 if power is negative')
    if p1 == 0:
        raise ValueError('Power law distribution undefined at -1')

    y=random.random()

    a = ((high**p1 - low**p1)*y + low**p1)**(1/p1)

    return a

def redshift_w(z=0, cosmology=True, w_min=0.1, w_max=5):
    """
    A random value from a uniform distribution redshifted by z

    Args:
        z (float): Redshift. Defaults to 0
        cosmology (boolean): Whether to use cosmology. Default to True
        w_min (float): Minimum pulse width [ms]. Defaults to 0.1
        w_max (float): Maximum pulse width [ms]. Defaults to 5
    Returns:
        w_int (float): A pulse width drawn from a uniform distribution, and
                       redshifted if cosmology was used.
    """
    w = random.uniform(w_min, w_max)
    if cosmology:
        w *= (1+z)
    return w
