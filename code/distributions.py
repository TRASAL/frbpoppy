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
