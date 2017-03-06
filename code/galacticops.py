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
    a_ngp = math.radians(12.9406333 * 15.)
    d_ngp = math.radians(27.1282500)
    l_ngp = math.radians(123.9320000)

    sd_ngp = math.sin(d_ngp)
    cd_ngp = math.cos(d_ngp)
    sb = math.sin(gb)
    cb = math.cos(gb)

    # Calculate right ascension
    y = cb*math.sin(l_ngp - gl)
    x = cd_ngp*sb - sd_ngp*cb*math.cos(l_ngp - gl)
    ra = math.atan2(y, x) + a_ngp
    ra = math.degrees(ra) % 360

    # Calculate declination
    dec = math.asin(sd_ngp*sb + cd_ngp*cb*math.cos(l_ngp - gl))
    dec = math.degrees(dec) % 360.
    if dec > 270:
        dec = - (360 - dec)
    return ra, dec


def calc_d_sun(x, y, z):
    """
    Calculate the distance from the source to the Sun

    Args:
        x (float): Galactic coordinate [kpc]
        y (float): Galactic coordinate [kpc]
        z (float): Galactic coordinate [kpc]
    Returns:
        d (float): Distance from source at the given galactic coordinates to
                   the Sun [kpc]
    """
    r_sun = 8.5  # kpc
    return math.sqrt(x**2 + (r_sun - y)**2 + z**2)


def ne2001_dist_to_dm(dist, gl, gb):
    """
    Convert position to a dispersion measure using NE2001 (compiled from
    fortran)

    Args:
        dist (float): Distance to source [kpc]. Distance will be cut at 100 kpc,
                      as NE2001 can not cope with larger distances. This value
                      should be more than enough to clear the Milky Way.
        gl (float): Galactic longitude [fractional degrees]
        gb (float): Galactic latitude [fractional degrees]
    Returns:
        dm (float): Dispersion measure [pc*cm^-3]
    """

    # NE2001 gives errors if distance input is too large! 100 kpc ought to be
    # enough to clear the galaxy.
    if dist > 100:
        dist = 100

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

def ne2001_get_smtau(dist, gl, gb):
    """
    Use the NE2001 model to calculate scattering measure. Calculations based on
    work presented in Cordes & Lazio (1991, DOI: 10.1086/170261)

    Args:
        dist (float): Distance to source [kpc]. Distance will be cut at 100 kpc,
                      as NE2001 can not cope with larger distances. Therefore
                      the calculated scattering will only be that from the
                      Milky Way.
        gl (float): Galactic longitude [fractional degrees]
        gb (float): Galactic latitude [fractional degrees]
    Returns:
        sm (float): Scattering measure
        smtau (float): Scattering measure, but unsure why different to sm
    """

    # NE2001 gives errors if distance input is too large! 100 kpc ought to be
    # enough to clear the galaxy.
    if dist > 100:
        dist = 100

    dist = C.c_float(dist)
    # Note the galactic coordinates need to be given in radians
    gl = C.c_float(math.radians(gl))
    gb = C.c_float(math.radians(gb))

    ndir = C.c_int(-1)
    sm = C.c_float(0.)
    smtau = C.c_float(0.)

    inpath = C.create_string_buffer(dm_mods.encode())
    linpath = C.c_int(len(dm_mods))

    ne2001lib.dmdsm_(C.byref(gl),
                     C.byref(gb),
                     C.byref(ndir),
                     C.byref(C.c_float(0.0)),
                     C.byref(dist),
                     C.byref(C.create_string_buffer(' '.encode())),
                     C.byref(sm),
                     C.byref(smtau),
                     C.byref(C.c_float(0.0)),
                     C.byref(C.c_float(0.0)),
                     C.byref(inpath),
                     C.byref(linpath)
                     )

    return sm.value, smtau.value


def ne2001_scint_time_bw(dist, gl, gb, freq):
    """
    Use the NE2001 model to get the diffractive scintillation timescale

    Args:
        dist (float): Distance to source [kpc]. Distance will be cut at 100 kpc,
                      as NE2001 can not cope with larger distances. Therefore
                      the calculated scintillation timescale will only be that
                      from the Milky Way.
        gl (float): Galactic longitude [fractional degrees]
        gb (float): Galactic latitude [fractional degrees]
        freq (float): Observing frequency [MHz]
    Returns:
        scint_time (float): Diffractive scintillation timescale [Hz]
        scint_bw (float): Scintillation bandwidth [Hz]
    """

    sm, smtau = ne2001_get_smtau(dist, gl, gb)

    if smtau <= 0.:
        scint_time = None
    else:
        # Eq. 46 of Cordes & Lazio 1991, ApJ, 376, 123 uses coefficient 3.3
        # instead of 2.3. They do this in the code and mention it explicitly,
        # so I trust it! <- From psrpoppy
        scint_time = 3.3 * (freq/1e3)**1.2 * smtau**(-0.6)

    if sm <= 0.:
        scint_bw = None
    else:
        # (eq. 48)
        scint_bw = 223. * (freq/1e3)**4.4 * sm**(-1.2) / dist

    return scint_time, scint_bw


def scatter_bhat(dm, scindex=-3.86, freq=1400.0):
    """
    Calculate scattering timescale according to Bhat et al. (2004,
    DOI:10.1086/382680) and to simluate the scatter around this relationship,
    draw from a Gaussian around this value.

    Args:
        dm (float): Dispersion measure [pc*cm^-3]
        scindex (float): Scattering index. Defaults to -3.86
        freq (float): Frequency at which to evaluate scattering time [MHz].
                      Defaults to 1400 MHz
    Returns:
        t_scat (float): Scattering timescale [ms]
    """

    log_t = -6.46 + 0.154*math.log10(dm) + 1.07*math.log10(dm)**2
    log_t += scindex*math.log10(freq/1e3)

    # Width of Gaussian distribution based on values given Lorimer et al. (2008)
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


def z_to_d(z, H_0=69.6, W_m=0.286, W_v=0.714):
    """
    Convert redshift to a luminosity distance. Based on James Schombert's python
    implementation of Edward L. Wright's cosmology calculator.

    Args:
        z (float): Redshift
        H_0 (float, optional): Hubble parameter. Defaults to 69.6
        W_m (float, optional): Omega matter. Defaults to 0.286
        W_k (float, optional): Omega vacuum. Defaults to 0.714
    Returns:
        d (float): Luminosity distance from Earth to the source [kpc]
    """

    # Initialize constants
    W_r = 0.4165/(H_0*H_0)  # Omega radiation
    W_k = 1.0 - W_m - W_r - W_v  # Omega curvature
    c = 299792.458  # Velocity of light [km/sec]
    dcmr = 0.
    az = 1/(1+z)

    # Calculate comoving distance
    n = 1000

    for i in range(n):
        a = az+(1-az)*(i+0.5)/n
        adot = math.sqrt(W_k + W_m/a + W_r/(a*a) + W_v*a*a)
        dcmr += 1/(a*adot)

    dcmr = (1.-az)*dcmr/n
    dc_mpc = (c/H_0)*dcmr  # Comoving distance [Mpc]

    # Calculate tangential comoving distance
    ratio = 1.
    x = math.sqrt(abs(W_k))*dcmr

    if x > 0.1:
        if W_k > 0:
            ratio = 0.5*(math.exp(x)-math.exp(-x))/x
        else:
            ratio = math.sin(x)/x
    else:
        y = x*x
        if W_k < 0:
            y = -y
        ratio = 1. + y/6. + y*y/120.

    dcmt = ratio*dcmr
    da = az*dcmt
    dl = da/(az*az)
    dl_mpc = (c/H_0)*dl # Luminosity distance [Mpc]

    return dl_mpc*1000


def d_to_z(d):
    """
    Convert distance in kpc to z. Holds for z <= 2

    Formulas from 'An Introduction to Modern Astrophysics (2nd Edition)' by
    Bradley W. Carroll, Dale A. Ostlie.

    Args:
        d (float): Distance from Earth to a source [kpc]
    Returns:
        z (float): Redshift
    """
    # Rewrite eq. 27.7 for z
    c = 299792.458  # Speed of light [km/s]
    H = 71.9  # Hubble's constant [km/s/Mpc] (Hubble Space Telescope)
    d /= 1e3  # Convert distance into units of Mpc
    dhc = d*H/c
    det = math.sqrt(1 - dhc**2)
    z = -(det + dhc - 1)/(dhc - 1)
    return z


def ioka_dm_igm(z):
    """
    Calculate the contribution of the intergalactic electron density to the
    dispersion measure, following Ioka (2003) and Inoue (2004)

    Args:
        z (float): Redshift of source
    Returns:
        dm_igm (float): Dispersion measure of the intergalactic medium [pc/cm^3]
    """
    return 1200 * z
