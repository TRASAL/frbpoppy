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
    y = random.random()
    p1 = power + 1
    return ((high**p1 - low**p1)*y + low**p1)**(1/p1)
