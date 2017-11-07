"""
Series of galactic operations (doesn't that sound cool?!), as in converting
coordinates, calculating DM etc.
"""

import ctypes as C
import csv
import math
import os
import random

# Import fortran libraries
mods = os.path.join(os.path.dirname(__file__), '../data/models/')
uni_mods = os.path.join(mods,'universe/')
dm_mods = os.path.join(mods, 'dm/')
loc = os.path.join(dm_mods, 'libne2001.so')
ne2001lib = C.CDLL(loc)
ne2001lib.dm_.restype = C.c_float


def lb_to_xyz(gl, gb, dist):
    """
    Convert galactic coordinates to galactic XYZ

    Args:
        l (float): Galactic longitude [fractional degrees]
        b (float): Galactic latitude [fractional degrees]
        dist (float): Distance to source [Gpc]

    Returns:
        gx, gy, gz: Galactic XYZ [Gpc]
    """
    rsun = 8.5e-6  # Gpc

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
    Bradley W. Carroll, Dale A. Ostlie (Eq. 24.19 onwards).

    NOTE: This function is not as accurate as the astropy conversion, nor as
    the Javascript calculators found online. However, as using astropy was
    prohibitively slow while running over large populations, frbpoppy uses this
    function. While this function is not as accurate, the under/over
    estimations of the coordinates are equally distributed meaning the errors
    cancel each other in the limit of large populations.

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
        dec = -(360 - dec)

    return ra, dec


def radec_to_lb(ra, dec):
    """
    Convert from ra, dec to galactic coordinates

    Formulas from 'An Introduction to Modern Astrophysics (2nd Edition)' by
    Bradley W. Carroll, Dale A. Ostlie (Eq. 24.16 onwards).

    NOTE: This function is not as accurate as the astropy conversion, nor as
    the Javascript calculators found online. However, as using astropy was
    prohibitively slow while running over large populations, frbpoppy uses this
    function. While this function is not as accurate, the under/over
    estimations of the coordinates are equally distributed meaning the errors
    cancel each other in the limit of large populations.

    Args:
        ra (string): Right ascension given in the form '19:06:53'
        dec (string): Declination given in the form '-40:37:14'

    Returns:
        gl, gb (float): Galactic longitude and latitude [fractional degrees]
    """

    # Convert to fractional degrees
    # Inspired by Joe Filippazzo calculator
    rh, rm, rs = [float(r) for r in ra.split(':')]
    ra = rh*15 + rm/4 + rs/240
    dd, dm, ds = [float(d) for d in dec.split(':')]
    if dd < 0:
        sign = -1
    else:
        sign = 1
    dec = dd + sign*dm/60 + sign*ds/3600

    a = math.radians(ra)
    d = math.radians(dec)

    # Coordinates of the galactic north pole (J2000)
    a_ngp = math.radians(12.9406333 * 15.)
    d_ngp = math.radians(27.1282500)
    l_ngp = math.radians(123.9320000)

    sd_ngp = math.sin(d_ngp)
    cd_ngp = math.cos(d_ngp)
    sd = math.sin(d)
    cd = math.cos(d)

    # Calculate galactic longitude
    y = cd*math.sin(a - a_ngp)
    x = cd_ngp*sd - sd_ngp*cd*math.cos(a - a_ngp)
    gl = - math.atan2(y, x) + l_ngp
    gl = math.degrees(gl) % 360
    # Shift so in range -180 to 180
    if gl > 180:
        gl = -(360 - gl)

    # Calculate galactic latitude
    gb = math.asin(sd_ngp*sd + cd_ngp*cd*math.cos(a - a_ngp))
    gb = math.degrees(gb) % 360.
    if gb > 270:
        gb = -(360 - gb)

    return gl, gb


def ergspers_to_watts(e):
    """Quick converstion from luminosity given in ergs/s to Watts"""
    return e*1e-7


def ne2001_dist_to_dm(dist, gl, gb):
    """
    Convert position to a dispersion measure using NE2001 (compiled from
    fortran)

    Args:
        dist (float): Distance to source [Gpc]. Distance will be cut at 100kpc,
                      as NE2001 can not cope with larger distances. This value
                      should be more than enough to clear the Milky Way.
        gl (float): Galactic longitude [fractional degrees]
        gb (float): Galactic latitude [fractional degrees]
    Returns:
        dm (float): Dispersion measure [pc*cm^-3]
    """
    dist *= 1e6  # Convert from Gpc to kpc

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
        dist (float): Distance to source [Gpc]. Distance will be cut at 100 kpc,
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
    dist *= 1e6  # Convert from Gpc to kpc

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


def scatter_bhat(dm, offset=-6.46, scindex=-3.86, freq=1400.0):
    """
    Calculate scattering timescale (values default to those from Bhat et al.
    (2004, DOI:10.1086/382680) and to simluate the scatter around this
    relationship, draw from a Gaussian around this value.

    Args:
        dm (float): Dispersion measure [pc*cm^-3]
        offset (float): Offset of scattering relationship. Defaults to -6.46
        scindex (float): Scattering index. Defaults to -3.86
        freq (float): Frequency at which to evaluate scattering time [MHz].
                      Defaults to 1400 MHz
    Returns:
        t_scat (float): Scattering timescale [ms]
    """

    log_t = offset + 0.154*math.log10(dm) + 1.07*math.log10(dm)**2
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


def z_to_v(z, H_0=69.6, W_m=0.286, W_v=0.714):
    """
    Convert redshift to a comoving volume. Based on James Schombert's python
    implementation of Edward L. Wright's cosmology calculator.

    Args:
        z (float): Redshift
        H_0 (float, optional): Hubble parameter. Defaults to 69.6
        W_m (float, optional): Omega matter. Defaults to 0.286
        W_k (float, optional): Omega vacuum. Defaults to 0.714
    Returns:
        v_gpc (float): Comoving volume from Earth to the source [Gpc^3]
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

    # Not necessary, but handy to have
    # Calculate luminosity distance
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
    dl_mpc = (c/H_0)*dl  # Luminosity distance [Mpc]

    # Calculate comoving volume
    ratio = 1.00
    x = math.sqrt(abs(W_k))*dcmr
    if x > 0.1:
        if W_k > 0:
            ratio = (0.125*(math.exp(2.*x)-math.exp(-2.*x))-x/2.)/(x*x*x/3.)
        else:
            ratio = (x/2. - math.sin(2.*x)/4.)/(x*x*x/3.)
    else:
        y = x*x
        if W_k < 0:
            y = -y
        ratio = 1. + y/5. + (2./105.)*y*y

    v_cm = ratio*dcmr*dcmr*dcmr/3.
    v_gpc = 4.*math.pi*((0.001*c/H_0)**3)*v_cm  # Comoving volume

    return v_gpc


def dist_lookup(cosmology=True, H_0=69.6, W_m=0.286, W_v=0.714, z_max=8.0):
    """
    Create a list of tuples to lookup the corresponding redshift for a comoving
    distance [Gpc]. Uses formulas from Hoggs et al. (1999) for the cosmological
    calculations, assuming a flat universe. To avoid long calculation times,
    it will check if a previous run with the same parameters has been done,
    which it will then load it. If not, it will calculate a new table, and save
    the table for later runs.

    Args:
        cosmology (boolean, optional): Whether to use cosmology or not.
        H_0 (float, optional): Hubble parameter. Defaults to 69.6
        W_m (float, optional): Omega matter. Defaults to 0.286
        W_k (float, optional): Omega vacuum. Defaults to 0.714
        z_max (float, optional): Maximum redshift. Defaults to 8.0
    Returns:
        None: If cosmology is not enabled
        ds (list): Comoving distances [Gpc]
        zs (list): Corresponding redshifts
    """
    # Initializing
    c = 299792.458  # Velocity of light [km/sec]
    ds = []
    zs = []

    if not cosmology:
        return None, None

    def cvt(value):
        """Convert a value to a string without a period"""
        return str(value).replace('.','d')

    # Filename
    paras = ['h0', cvt(H_0),
             'wm', cvt(W_m),
             'wv', cvt(W_v),
             'zmax', cvt(z_max)]
    f = '-'.join(paras) + '.csv'

    # Check whether frbpoppy can avoid creating new tables
    if os.path.isfile(uni_mods + f):
        with open(uni_mods + f, 'r') as f:
            data = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
            for r in data:
                ds.append(r[0])
                zs.append(r[1])

    # If unavoidable, proceed with integration...
    else:

        # Knowingly put imports here
        from scipy.integrate import quad as integrate
        import numpy as np

        # Initialize parameters
        step = 0.001
        W_k = 1.0 - W_m - W_v  # Omega curvature

        if W_k != 0.0:
            print('Careful - Your cosmological parameters do not sum to 1.0')

        # Numerically integrate the following function
        def d_c(x):
            """Comoving distance (Hogg et al, 1999)"""
            return 1/math.sqrt((W_m*(1+x)**3 + W_k*(1+x)**2 + W_v))

        for z in np.arange(0,z_max+step,step):
            d = c/H_0*integrate(d_c, 0, z)[0]
            d /= 1e3  # Covert from Mpc to Gpc
            ds.append(d)
            zs.append(z)

        # Save distance and redshift
        with open(uni_mods + f, 'w+') as df:
            d = '\n'.join('{},{}'.format(ds[i],zs[i]) for i in range(len(ds)))
            df.write(d)

    return ds, zs


def interpolate_z(d, ds, zs, H_0=69.6):
    """
    Interpolate between two comoving distances.

    Obtain an approximate redshift, unless dzs is none (i.e. cosmology has
    been set to False in the function dist_lookup), in which case use an
    approximation that works to z<2.

    Args:
        d (float): Comoving distance [Gpc]
        ds (list): Comoving distances [Gpc]
        zs (list): Corresponding redshifts
        H_0 (float, optional): The value of the Hubble expansion, only required
            when using a non-cosmology approximation. Defaults to 69.6
    Returns:
        z (float): Redshift

    """
    z = 0

    if not ds:
        # Convert distance in Gpc to z. Holds for z <= 2
        # Formulas from 'An Introduction to Modern Astrophysics (2nd Edition)'
        # by Bradley W. Carroll, Dale A. Ostlie.
        c = 299792.458  # Velocity of light [km/sec]
        d *= 1e3
        dhc = d*H_0/c
        det = math.sqrt(1 - dhc**2)
        z = -(det + dhc - 1)/(dhc - 1)
        return z

    # Interpolate between values
    for i,e in enumerate(ds):

        if d < e:
            i2 = i
            i1 = i-1

            # Return lowest z if before first value
            if i1 < 0:
                return zs[0]

            slope = (zs[i2] - zs[i1]) / (ds[i2] - ds[i1])
            z = zs[i1] + slope*(d - ds[i1])

            return z

    print('Gone over your maximum redshift', d)
    # Return highest z if beyond the largest redshift
    return zs[-1]


def ioka_dm_igm(z, slope=1200):
    """
    Calculate the contribution of the intergalactic electron density to the
    dispersion measure, following Ioka (2003) and Inoue (2004)

    Args:
        z (float): Redshift of source
        slope (int, optional): Slope of relationship
    Returns:
        dm_igm (float): Dispersion measure of intergalactic medium [pc/cm^3]
    """
    return random.gauss(slope*z, 0.2*slope*z)
