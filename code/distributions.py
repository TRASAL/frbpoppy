import math
import numpy as np
import numpy.fft as fft
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
        cosmology (boolean): Whether to use cosmology
        w_min (float): Minimum pulse width [ms]
        w_max (float): Maximum pulse width [ms]
    Returns:
        w_int (float): A pulse width drawn from a uniform distribution, and
                       redshifted if cosmology was used.
    """
    w = random.uniform(w_min, w_max)
    if cosmology:
        w *= (1+z)
    return w


def pink_noise():
    """
    Simluate burst times using pink noise

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
    Following Oppermann & Pen (2017), simulate repeat times

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

    return ts[:-1]