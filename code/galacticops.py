"""
Series of galactic operations (doesn't that sound cool?!), as in converting
coordinates, calculating DM etc.
"""

import math


def lb_to_xyz(gl, gb, dist):
    """
    Convert galactic coordinates to galactic XYZ

    Args:
        dist (float): Distance to source in kpc
    """

    if gl is None or gb is None:
        raise IOError('Galactic longitude and/or latitude not set')

    rsun = 8.5  # kpc

    L = math.radians(gl)
    B = math.radians(gb)

    gx = dist * math.cos(B) * math.sin(L)
    gy = rsun - dist * math.cos(B) * math.cos(L)
    gz = dist * math.sin(B)

    return (gx, gy, gz)


def lb_to_radec(l, b):
    """
    Convert galactic coordinates to RA, Dec

    Formulas from 'An Introduction to Modern Astrophysics (2nd Edition)' by
    Bradley W. Carroll, Dale A. Ostlie.

    Args:
        l (float): Galactic longitude in fractional degrees
        b (float): Galactic latitude in fractional degrees

    Returns:
        ra, dec (float): Right ascension and declination in fractional degrees
    """

    gl = math.radians(l)
    gb = math.radians(b)

    # Coordinates of the galactic north pole (J2000)
    a_ngp = math.radians(12.9406333 * 15. + 180.)
    d_ngp = math.radians(27.1282500)
    l_ngp = math.radians(123.9320000)

    sd_ngp = math.sin(d_ngp)
    cd_ngp = math.cos(d_ngp)
    sb = math.sin(gb)
    cb = math.cos(gb)

    # Calculate right ascension
    y = cb*math.sin(l_ngp - gl)
    x = cd_ngp*sb - sd_ngp*cb*math.cos(l_ngp - gl)
    ra = math.atan(y/x) + a_ngp
    ra = math.degrees(ra)

    # Ensure value is in right quadrant
    if y > 0 and x > 0:
        ra = ra % 90
    elif y > 0 and x < 0:
        ra = 90 + (ra % 90)
    elif y < 0 and x < 0:
        ra = 180 + (ra % 90)
    elif y < 0 and x > 0:
        ra = 270 + (ra % 90)

    # Calculate declination
    dec = math.asin(sd_ngp*sb + cd_ngp*cb*math.cos(l_ngp - gl))
    dec = math.degrees(dec) % 360.

    return ra, dec
