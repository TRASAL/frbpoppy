"""
Series of galactic operations (doesn't that sound cool?!).

...as in converting coordinates, calculating DM etc.
"""
from datetime import timedelta
import ctypes as C
import math
import numpy as np
import os
import pandas as pd
import random

from frbpoppy.paths import paths

# Import fortran libraries
uni_mods = os.path.join(paths.models(), 'universe/')
dm_mods = os.path.join(paths.models(), 'ne2001/')
loc = os.path.join(dm_mods, 'libne2001.so')
ne2001lib = C.CDLL(loc)
ne2001lib.dm_.restype = C.c_float


def frac_deg(ra, dec):
    """Convert coordinates expressed in hh:mm:ss to fractional degrees."""
    # Inspired by Joe Filippazzo calculator
    rh, rm, rs = [float(r) for r in ra.split(':')]
    ra = rh*15 + rm/4 + rs/240
    dd, dm, ds = [float(d) for d in dec.split(':')]
    if dd < 0:
        sign = -1
    else:
        sign = 1
    dec = dd + sign*dm/60 + sign*ds/3600
    return ra, dec


def lb_to_xyz(gl, gb, dist):
    """
    Convert galactic coordinates to galactic XYZ.

    Args:
        l (float): Galactic longitude [fractional degrees]
        b (float): Galactic latitude [fractional degrees]
        dist (float): Distance to source [Gpc]

    Returns:
        gx, gy, gz: Galactic XYZ [Gpc]

    """
    rsun = 8.5e-6  # Gpc

    L = np.radians(gl)
    B = np.radians(gb)

    gx = dist * np.cos(B) * np.sin(L)
    gy = rsun - dist * np.cos(B) * np.cos(L)
    gz = dist * np.sin(B)

    return gx, gy, gz


def lb_to_radec(l, b):
    """
    Convert galactic coordinates to RA, Dec.

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
    gl = np.radians(l)
    gb = np.radians(b)

    # Coordinates of the galactic north pole (J2000)
    a_ngp = np.radians(12.9406333 * 15.)
    d_ngp = np.radians(27.1282500)
    l_ngp = np.radians(123.9320000)

    sd_ngp = np.sin(d_ngp)
    cd_ngp = np.cos(d_ngp)
    sb = np.sin(gb)
    cb = np.cos(gb)

    # Calculate right ascension
    y = cb*np.sin(l_ngp - gl)
    x = cd_ngp*sb - sd_ngp*cb*np.cos(l_ngp - gl)
    ra = np.arctan2(y, x) + a_ngp
    ra = np.degrees(ra) % 360

    # Calculate declination
    dec = np.arcsin(sd_ngp*sb + cd_ngp*cb*np.cos(l_ngp - gl))
    dec = np.degrees(dec) % 360.
    dec[dec > 270] = -(360 - dec[dec > 270])

    return ra, dec


def radec_to_lb(ra, dec, frac=False):
    """
    Convert from ra, dec to galactic coordinates.

    Formulas from 'An Introduction to Modern Astrophysics (2nd Edition)' by
    Bradley W. Carroll, Dale A. Ostlie (Eq. 24.16 onwards).

    NOTE: This function is not as accurate as the astropy conversion, nor as
    the Javascript calculators found online. However, as using astropy was
    prohibitively slow while running over large populations, we use this
    function. While this function is not as accurate, the under/over
    estimations of the coordinates are equally distributed meaning the errors
    cancel each other in the limit of large populations.

    Args:
        ra (string): Right ascension given in the form '19:06:53'
        dec (string): Declination given in the form '-40:37:14'
        frac (bool): Denote whether coordinates are already fractional or not
    Returns:
        gl, gb (float): Galactic longitude and latitude [fractional degrees]

    """
    if not frac:
        ra, dec = frac_deg(ra, dec)

    a = np.radians(ra)
    d = np.radians(dec)

    # Coordinates of the galactic north pole (J2000)
    a_ngp = np.radians(12.9406333 * 15.)
    d_ngp = np.radians(27.1282500)
    l_ngp = np.radians(123.9320000)

    sd_ngp = np.sin(d_ngp)
    cd_ngp = np.cos(d_ngp)
    sd = np.sin(d)
    cd = np.cos(d)

    # Calculate galactic longitude
    y = cd*np.sin(a - a_ngp)
    x = cd_ngp*sd - sd_ngp*cd*np.cos(a - a_ngp)
    gl = - np.arctan2(y, x) + l_ngp
    gl = np.degrees(gl) % 360

    # Shift so in range -180 to 180
    if isinstance(gl, np.ndarray):
        gl[gl > 180] = -(360 - gl[gl > 180])
    else:
        if gl > 180:
            gl = -(360 - gl)

    # Calculate galactic latitude
    gb = np.arcsin(sd_ngp*sd + cd_ngp*cd*np.cos(a - a_ngp))
    gb = np.degrees(gb) % 360

    if isinstance(gb, np.ndarray):
        gb[gb > 270] = -(360 - gb[gb > 270])
    else:
        if gb > 270:
            gb = -(360 - gb)

    return gl, gb


def separation(ra_1, dec_1, ra_2, dec_2):
    """Separation between points on sky [degrees].

    Using a special case of the Vincenty formula for an ellipsoid with equal
    major and minor axes.

    See https://en.wikipedia.org/wiki/Great-circle_distance for more info.

    """
    # Convert to radians
    ra_1 = np.deg2rad(ra_1)
    dec_1 = np.deg2rad(dec_1)
    ra_2 = np.deg2rad(ra_2)
    dec_2 = np.deg2rad(dec_2)

    # Shortcuts
    sdr = np.sin(ra_2 - ra_1)
    cdr = np.cos(ra_2 - ra_1)
    cd1 = np.cos(dec_1)
    cd2 = np.cos(dec_2)
    sd1 = np.sin(dec_1)
    sd2 = np.sin(dec_2)

    # Calculation
    upper = np.sqrt((cd2*sdr)**2 + (cd1*sd2 - sd1*cd2*cdr)**2)
    lower = sd1*sd2 + cd1*cd2*cdr
    sep = np.arctan2(upper, lower)

    return np.rad2deg(sep)


def ergspers_to_watts(e):
    """Quick converstion from luminosity given in ergs/s to Watts."""
    return e*1e-7


def ne2001_dist_to_dm(dist, gl, gb):
    """
    Convert position to a dispersion measure using NE2001.

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
    Use the NE2001 model to calculate scattering measure.

    Calculations based on work presented in Cordes & Lazio
    (1991, DOI: 10.1086/170261)

    Args:
        dist (array): Distance to source [kpc]. Distance will be cut at 100 kpc
                      as NE2001 can not cope with larger distances. Therefore
                      the calculated scattering will only be that from the
                      Milky Way.
        gl (array): Galactic longitude [fractional degrees]
        gb (array): Galactic latitude [fractional degrees]
    Returns:
        sm (array): Scattering measure
        smtau (array): Scattering measure, but unsure why different to sm

    """
    # NE2001 gives errors if distance input is too large! 100 kpc ought to be
    # enough to clear the galaxy.
    dist[dist > 100] = 100
    sms = np.ones_like(dist)
    smtaus = np.ones_like(dist)

    for i, d in enumerate(dist):

        disti = C.c_float(d)
        # Note the galactic coordinates need to be given in radians
        gli = C.c_float(math.radians(gl[i]))
        gbi = C.c_float(math.radians(gb[i]))

        ndir = C.c_int(-1)
        sm = C.c_float(0.)
        smtau = C.c_float(0.)

        inpath = C.create_string_buffer(dm_mods.encode())
        linpath = C.c_int(len(dm_mods))

        ne2001lib.dmdsm_(C.byref(gli),
                         C.byref(gbi),
                         C.byref(ndir),
                         C.byref(C.c_float(0.0)),
                         C.byref(disti),
                         C.byref(C.create_string_buffer(' '.encode())),
                         C.byref(sm),
                         C.byref(smtau),
                         C.byref(C.c_float(0.0)),
                         C.byref(C.c_float(0.0)),
                         C.byref(inpath),
                         C.byref(linpath)
                         )

        sms[i], smtaus[i] = sm.value, smtau.value

    return sms, smtaus


def ne2001_scint_time_bw(dist, gl, gb, freq):
    """
    Use the NE2001 model to get the diffractive scintillation timescale.

    Args:
        dist (array): Distance to source [Gpc]. Distance will be cut at 100 kpc
                      as NE2001 can not cope with larger distances. Therefore
                      the calculated scintillation timescale will only be that
                      from the Milky Way.
        gl (array): Galactic longitude [fractional degrees]
        gb (array): Galactic latitude [fractional degrees]
        freq (float): Observing frequency [MHz]
    Returns:
        scint_time (float): Diffractive scintillation timescale [Hz]
        scint_bw (float): Scintillation bandwidth [Hz]

    """
    dist *= 1e6  # Convert from Gpc to kpc

    sm, smtau = ne2001_get_smtau(dist, gl, gb)

    scint_time = np.ones_like(dist)
    scint_time[smtau <= 0.] = float('NaN')
    # Eq. 46 of Cordes & Lazio 1991, ApJ, 376, 123 uses coefficient 3.3
    # instead of 2.3. They do this in the code and mention it explicitly,
    # so I trust it! <- From psrpoppy
    scint_time[smtau > 0.] = 3.3 * (freq/1e3)**1.2 * smtau**(-0.6)

    scint_bw = np.ones_like(dist)
    scint_bw[sm <= 0.] = float('NaN')
    # (eq. 48)
    scint_bw[sm > 0.] = 223. * (freq/1e3)**4.4 * sm**(-1.2) / dist

    return scint_time, scint_bw


def scatter_bhat(dm, offset=-6.46, scindex=-3.86, freq=1400.0):
    """
    Calculate scattering timescale (values default to those from Bhat et al.
    (2004, DOI:10.1086/382680) and to simluate the scatter around this
    relationship, draw from a Gaussian around this value.

    Args:
        dm (array): Dispersion measure [pc*cm^-3]
        offset (float): Offset of scattering relationship. Defaults to -6.46
        scindex (float): Scattering index. Defaults to -3.86
        freq (float): Frequency at which to evaluate scattering time [MHz].
                      Defaults to 1400 MHz
    Returns:
        array: Scattering timescale [ms]

    """
    log_t = offset + 0.154*np.log10(dm) + 1.07*np.log10(dm)**2
    log_t += scindex*np.log10(freq/1e3)

    # Width of Gaussian distribution based on values given Lorimer et al (2008)
    t_scat = 10**np.random.normal(log_t, 0.8)

    return t_scat


def load_T_sky():
    """
    Read the Haslam sky temperature map into a list.

    ... from which temperatures can
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
                except ValueError:
                    pass
                str_idx += 5

    return t_sky_list


class Redshift:
    """Class for converting redshift to other distance measures."""

    def __init__(self, z, H_0=67.74, W_m=0.3089, W_v=0.6911):
        """
        Convert redshift to a various measures.

        Based on James Schombert's python implementation of Edward L. Wright's
        cosmology calculator.

        Args:
            z (array): Redshift
            self.H_0 (float, optional): Hubble parameter.
            self.W_m (float, optional): Omega matter.
            self.W_v (float, optional): Omega vacuum.

        Returns:
            array: One of the distance measures [Gpc], or comoving volume from
                Earth to the source [Gpc^3]

        """
        self.z = z
        self.H_0 = H_0
        self.W_m = W_m
        self.W_v = W_v

        # Initialize constants
        self.W_r = 0.4165/(self.H_0*self.H_0)  # Omega radiation
        self.W_k = 1.0 - self.W_m - self.W_r - self.W_v  # Omega curvature
        self.c = 299792.458  # Velocity of light [km/sec]
        self.dcmr = 0.
        self.az = 1/(1+self.z)

        # Distance measures
        self.dc_mpc = None
        self.dl_mpc = None

    def dist_co(self):
        """Calculate the corresponding comoving distance [Gpc]."""
        n = 1000

        for i in range(n):
            a = self.az+(1-self.az)*(i+0.5)/n
            s = sum([self.W_k, self.W_m/a, self.W_r/(a*a), self.W_v*a*a])
            adot = np.sqrt(s)
            self.dcmr += 1/(a*adot)

        self.dcmr = (1.-self.az)*self.dcmr/n

        self.dc_mpc = (self.c/self.H_0)*self.dcmr  # Comoving distance [Mpc]

        return self.dc_mpc*1e-3  # Convert to Gpc

    def dist_lum(self):
        """Calculate the corresponding luminosity distance [Gpc]."""
        if self.dc_mpc is None:
            self.dist_co()

        # Calculate luminosity distance
        ratio = np.ones_like(self.dcmr)
        x = np.sqrt(abs(self.W_k))*self.dcmr

        mask = (x > 0.1)

        if self.W_k > 0:
            ratio[mask] = 0.5*(np.exp(x[mask])-np.exp(-x[mask]))/x[mask]
        else:
            ratio[mask] = np.sin(x[mask])/x[mask]

        y = x*x
        if self.W_k < 0:
            y = -y
        ratio[~mask] = 1. + y[~mask]/6. + y[~mask]*y[~mask]/120.

        dcmt = ratio*self.dcmr
        da = self.az*dcmt
        dl = da/(self.az*self.az)
        self.dl_mpc = (self.c/self.H_0)*dl  # Luminosity distance [Mpc]

        return self.dl_mpc*1e-3  # Covert to Gpc

    def vol_co(self):
        """Calculate the corresponding comoving volume [Gpc^3]."""
        if self.dl_mpc is None:
            self.dist_lum()

        ratio = np.ones_like(self.dcmr)
        x = math.sqrt(abs(self.W_k))*self.dcmr

        mask = (x > 0.1)

        if self.W_k > 0:
            n = (0.125*(np.exp(2.*x[mask])-np.exp(-2.*x[mask]))-x[mask]/2.)
            ratio[mask] = n/(x[mask]**3/3)
        else:
            ratio[mask] = (x[mask]/2. - np.sin(2.*x[mask])/4.)/(x[mask]**3/3)

        y = x*x
        if self.W_k < 0:
            y = -y
        ratio[~mask] = 1. + y[~mask]/5. + (2./105.)*y[~mask]*y[~mask]

        v_cm = ratio*self.dcmr**3/3
        self.v_gpc = 4.*math.pi*((1e-3*self.c/self.H_0)**3)*v_cm

        return self.v_gpc


def z_to_d_approx(z, H_0=67.74):
    """
    Calculate distance in Gpc from a redshift.

    Only holds for z <= 2. Formulas from 'An Introduction to Modern
    Astrophysics (2nd Edition)' by Bradley W. Carroll, Dale A. Ostlie.
    (Eq. 27.7)

    Args:
        z (float): Redshift
        H_0 (float, optional): Hubble parameter. Defaults to 67.74
    Returns:
        dist (float): Associated distance [Gpc]
    """
    c = 299792.458  # Velocity of light [km/sec]
    zsq = (z+1)**2
    dist = c/H_0 * (zsq - 1)/(zsq + 1)
    dist /= 1e3  # Mpc -> Gpc
    return dist


def dist_to_z(dist, H_0=67.74):
    """
    Calculate redshift from a distance in Gpc.

    Only holds for z <= 2. Formulas from 'An Introduction to Modern
    Astrophysics (2nd Edition)' by Bradley W. Carroll, Dale A. Ostlie.
    (Eq. 27.7)

    Args:
        dist (float): Distance [Gpc].
        H_0 (float, optional): Hubble parameter. Defaults to 67.74
    Returns:
        z (float): Associated redshift
    """
    c = 299792.458  # Velocity of light [km/sec]
    dist *= 1e3  # Gpc -> Mpc
    dhc = dist*H_0/c
    det = math.sqrt(1 - dhc**2)
    z = -(det + dhc - 1)/(dhc - 1)
    return z


def datetime_to_julian(date):
    """Convert a datetime object into julian float.

    See https://aa.usno.navy.mil/faq/docs/JD_Formula.php for more info.

    Args:
        date (datetime-object): Date in question

    Returns:
        float: Julian calculated datetime.

    """
    # Add support for numpy arrays of datetime64
    if np.issubdtype(date.dtype, np.datetime64):
        date = pd.to_datetime(date)

    # Define terms
    y = date.year
    m = date.month
    d = date.day
    h = date.hour
    min = date.minute
    sec = date.second

    # Calculate julian day number
    jdn = 367*y - ((7*(y + ((m+9)/12).astype(int)))/4).astype(int)
    jdn += ((275*m)/9).astype(int) + d + 1721013.5
    # Add fractional day
    jd = jdn + h/24 + min/1440 + sec/86400

    # Convert to a numpy array
    if isinstance(jd, pd.Float64Index):
        jd = jd.values

    return jd


def datetime_to_gmst(date):
    """Calculate Greenwich Mean Sidereal Time.

    See https://aa.usno.navy.mil/faq/docs/GAST.php for more info.
    """
    jd = datetime_to_julian(date)
    return ((18.697374558 + 24.06570982441908*(jd - 2451545))*15) % 360


def random_date(start, end):
    """Generate a random datetime between two datetime objects."""
    delta = end - start
    int_delta = (delta.days * 24 * 60 * 60) + delta.seconds
    random_second = random.randrange(int_delta)
    return start + timedelta(seconds=random_second)


def coord_to_offset(xref, yref, x, y):
    """
    Convert point (x, y) to projected offset from reference (xref, yref).

    Makes use of a gnomonic projection: see both
        https://github.com/LSSTDESC/Coord/blob/master/coord/celestial.py
        http://mathworld.wolfram.com/GnomonicProjection.html

    Args:
        xref (array): Reference RA or Az [rad]
        yref (array): Reference Dec or Alt [rad]
        x (array): Target RA or Az [rad]
        y (array): Target Dec or Alt [rad]

    Returns:
        array, array: x and y offset [rad]

    """
    # Define convenience numbers
    sinxref = np.sin(xref)
    sinx = np.sin(x)
    cosxref = np.cos(xref)
    cosx = np.cos(x)

    sinyref = np.sin(yref)
    siny = np.sin(y)
    cosyref = np.cos(yref)
    cosy = np.cos(y)

    # Sine and cosine of shift in x
    cosdx = (cosxref * cosx) + (sinxref * sinx)
    sindx = (cosxref * sinx) - (sinxref * cosx)

    # Projection effect cosine
    cosc = sinyref * siny + cosyref * cosy * cosdx

    # Projected offsets
    dx = (cosy * sindx) / cosc
    dy = (cosyref * siny - sinyref * cosy * cosdx) / cosc

    if cosc < 0:
        dx, dy = np.nan, np.nan

    return dx, dy


def hadec_to_azalt(ha, dec, lat):
    """
    Convert hour angle and declination to azimuth and altitude.

    Args:
        ha (array): Hour angle [rad]
        dec (array): Declination [rad]
        lat (float): Latitude [rad]

    Returns:
        array, array: az, alt [rad]

    """
    # Ha and dec should be same type
    assert type(ha) == type(dec)

    # Altitude
    sinalt = np.sin(dec) * np.sin(lat) + np.cos(dec) * np.cos(lat) * np.cos(ha)
    alt = np.arcsin(sinalt)

    # azimuth (note this uses altitude)
    cosaz = (np.sin(dec)-np.sin(alt)*np.sin(lat)) / (np.cos(alt)*np.cos(lat))

    convert_to_float = False
    if isinstance(cosaz, float):
        cosaz = np.array([cosaz])
        convert_to_float = True

    # Numerical instability can cause cosaz > 1
    cosaz[cosaz > 1] = 1
    cosaz[cosaz < -1] = -1
    az = np.arccos(cosaz)

    # Sign of azimuth is lost, but can be recovered using the input hour angle
    mask = np.sin(ha) > 0
    az[mask] = 2*np.pi - az[mask]

    # Convert back to float if input was float
    if convert_to_float:
        az = float(az)

    return az, alt


def in_region(ra, dec, gl, gb,
              ra_min=0, ra_max=360,
              dec_min=-90, dec_max=90,
              gl_min=-180, gl_max=180,
              gb_min=-90, gb_max=90):
    """
    Check if the given frbs are within the survey region.

    Args:
        ra, dec, gl, gb (float): Coordinates to check whether in region

    Returns:
        array: Boolean mask denoting whether frbs are within survey region

    """
    # Create mask with False
    mask = np.ones_like(ra, dtype=bool)

    # Ensure in correct format
    gl[gl > 180.] -= 360.

    # Create region masks
    gl_limits = (gl > gl_max) | (gl < gl_min)
    gb_limits = (gb > gb_max) | (gb < gb_min)
    ra_limits = (ra > ra_max) | (ra < ra_min)
    dec_limits = (dec > dec_max) | (dec < dec_min)

    mask[gl_limits] = False
    mask[gb_limits] = False
    mask[ra_limits] = False
    mask[dec_limits] = False

    return mask


def calc_sky_radius(area):
    """Determine the radius along the sky of an area in sq. degrees."""
    # Check whether the full sky
    if np.allclose(area, 4*np.pi*(180/np.pi)**2):
        return 180
    else:
        cos_r = (1 - (area*np.pi)/(2*180**2))
        return np.rad2deg(np.arccos(cos_r))


def calc_sky_area(radius):
    """Determine the area in sq. degree of a radius along the sky."""
    return (1 - np.cos(np.deg2rad(radius)))*(2*180**2)/np.pi
