"""
Series of galactic operations (doesn't that sound cool?!), as in converting
coordinates, calculating DM etc.
"""

import ctypes as C
import math
import os
import random

# Import fortran libraries
dm_mods = os.path.join(os.path.dirname(__file__), '../data/models/dm/')
loc = os.path.join(dm_mods, 'libne2001.so')
ne2001lib = C.CDLL(loc)
ne2001lib.dm_.restype = C.c_float


def lb_to_xyz(gl, gb, dist):
    """
    Convert galactic coordinates to galactic XYZ

    Args:
        dist (float): Distance to source [kpc]
    """

    if gl is None or gb is None:
        raise IOError('Galactic longitude and/or latitude not set')

    rsun = 8.5  # kpc

    L = math.radians(gl)
    B = math.radians(gb)

    gx = dist * math.cos(B) * math.sin(L)
    gy = rsun - dist * math.cos(B) * math.cos(L)
    gz = dist * math.sin(B)

    return gx, gy, gz


def lb_to_radec(l, b):
    """
    Convert galactic coordinates to RA, Dec

    Formulas from 'An Introduction to Modern Astrophysics (2nd Edition)' by
    Bradley W. Carroll, Dale A. Ostlie.

    Args:
        l (float): Galactic longitude [fractional degrees]
        b (float): Galactic latitude [fractional degrees]

    Returns:
        ra, dec (float): Right ascension and declination [fractional degrees]
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


def calc_d_sun(x, y, z):
    """Calculate the distance from the source to the Sun (all units kpc)"""
    r_sun = 8.5  # kpc
    return math.sqrt(x**2 + (r_sun - y)**2 + z**2)


def ne2001_dist_to_dm(dist, gl, gb):
    """
    Convert position to a dispersion measure using NE2001 (compiled from
    fortran)

    Args:
        dist (float): Distance to source [kpc]
        gl (float): Galactic longitude [fractional degrees]
        gb (float): Galactic latitude [fractional degrees]
    Returns:
        dm (float): Dispersion measure [pc*cm^-3]
    """

    dist = C.c_float(dist)
    gl = C.c_float(gl)
    gb = C.c_float(gb)
    inpath = C.create_string_buffer(dm_mods.encode())
    linpath = C.c_int(len(dm_mods))

    dm = ne2001lib.dm_(C.byref(dist),
                       C.byref(gl),
                       C.byref(gb),
                       C.byref(C.c_int(4)),
                       C.byref(C.c_float(0.0)),
                       C.byref(inpath),
                       C.byref(linpath)
                       )

    return dm


def scatter_bhat(dm, scindex=-3.86, freq=1400.0):
    """
    Calculate scattering timescale according to Bhat et al. (2004,
    DOI:10.1086/382680) and to simluate the scatter around this relationship,
    draw from a Gaussian around this value.

    Args:
        dm (float): Dispersion measure [pc*cm^-3]
        freq (float): Frequency at which to evaluate scattering time [MHz]
    Returns:
        t_scat (float): Scattering timescale [ms]
    """

    log_t = -6.46 + 0.154*math.log10(dm) + 1.07*math.log10(dm)**2
    log_t += scindex*math.log10(freq/1e3)

    # Width of Gaussian distribution based on values given in psrpoppy
    t_scat = 10**random.gauss(log_t, 0.8)

    return t_scat


def load_T_sky():
    """
    Read the Haslam sky temperature map into a list from which temperatures can
    be retrieved. The temperature sky map is given in the weird units of
    HealPix, and despite looking up info on this coordinate system, I don't
    have the foggiest idea of how to transform these to galactic coordinates. I
    have therefore directly copied the following code from psrpoppy in the
    assumption Sam Bates managed to figure it out.

    Returns:
        t_sky_list (list): List of sky temperatures in HealPix? coordinates?
    """

    model = os.path.join(os.path.dirname(__file__), '../data/models/tsky/')
    path = os.path.join(model, 'haslam_2014.dat')

    t_sky_list = []
    with open(path) as f:
        for line in f:
            str_idx = 0
            while str_idx < len(line):
                # each temperature occupies space of 5 chars
                temp_string = line[str_idx:str_idx+5]
                try:
                    t_sky_list.append(float(temp_string))
                except:
                    pass
                str_idx += 5

    return t_sky_list
