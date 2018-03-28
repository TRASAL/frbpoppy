"""Define distributions from which to get random numbers."""
import numpy as np
import random


def powerlaw(low, high, power):
    """
    Return random variable distributed according to power law.

    The power law distribution power is simply the power, not including a minus
    sign (P scales with x^n with n the power). A flat powerlaw can therefore
    be created by taking setting power to zero.

    Args:
        low (float): Lower limit of distribution
        high (float): Higher limit of distribution
        power (float): Power of power law distribution
    Returns:
        x (float): Random variable picked from power law distribution
    """
    if low > high:
        low, high = high, low

    # p1 = power
    if (low == 0. or high == 0.) and power < 0:
        raise ValueError('Power law not defined at 0 if power is negative')
    # Not completely kosher, but hey...
    if power == 0:
        power = 1e-15

    y = random.random()

    a = ((high**power - low**power)*y + low**power)**(1/power)

    return a


def pink_noise():
    """
    Simluate burst times using pink noise.

    Returns:
        ts (list): A list of burst times

    """
    # Assume FRBs can repeat at max once per 20s, and that
    # would be observable for a maximum of 12h
    length = round(86400 / 2 / 20)

    # Create gaussion noise
    signal = np.random.normal(0, 1, size=length)
    time = [i for i in range(length)]

    # Fast Fourier Transform it
    freq = np.fft.rfftfreq(length)
    four_trans = np.fft.rfft(signal)

    # Turn into pink noise
    pn = [e[0]*e[1]**-1 for e in zip(four_trans[1:], freq[1:])]

    # Revert to time domain
    ns = np.fft.irfft(pn)

    # If above limit, then frb
    limit = max(ns) - 2
    ts = [e[1] for e in zip(ns, time[:-2]) if e[0] > limit]

    # Revert to normal timescale
    ts *= 20

    return ts


def oppermann_pen():
    """
    Following Oppermann & Pen (2017), simulate repeat times.

    Returns:
        ts (list): List of burst times

    """
    r = 5.7
    k = 0.34

    ts = []
    t_tot = 0.5  # Assuming a maximum of 12 hours on one spot
    t_sum = 0.0

    # Get time of bursts
    while t_sum < t_tot:
        t = r*np.random.weibull(k)
        t_sum += t
        ts.append(t_sum)

    # Convert to seconds
    ts = [t*86400 for t in ts[:-1]]

    return ts


def z_from_sfr(z_max=2.5):
    """
    Return a random redshift for sources following the Star Formation Rate.

    Follows Madau & Dickinson (2014), eq. 15. For more info see
    https://arxiv.org/pdf/1403.0007.pdf

    """
    def sfr(z):
        return (1+z)**2.7/(1+((1+z)/2.9)**5.6)

    z = None

    while not z:
        x = random.random()*z_max
        y = random.random()*9.0
        if y <= sfr(x):
            z = x

    return z
