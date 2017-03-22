import random


def powerlaw(low, high, power, y=random.random()):
    """
    Return random variable distributed according to power law

    Args:
        low (float): Lower limit of distribution
        high (float): Higher limit of distribution
        power (float): Power of power law distribution
        y (float): Random seed variable
    Returns:
        x (float): Random variable picked from power law distribution
    """
    if low > high:
        low, high = high, low

    p1 = power + 1
    if (low == 0. or high == 0.) and p1 < 0:
        raise ValueError('Power laws are not defined at 0 if power is negative')

    return ((high**p1 - low**p1)*y + low**p1)**(1/p1)
