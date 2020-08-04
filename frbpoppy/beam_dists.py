"""Functions for calculating the intensity of points in a beam."""
import numpy as np
from scipy.special import j1

import frbpoppy.galacticops as go
from frbpoppy.paths import paths


def get_beam_props(model, fwhm):
    """Get beam properties.

    Args:
        model (str): Which model to use.
        fwhm (float): FWHM [frac. deg].

    Returns:
        beam_size, pixel_scale, beam_array

    """
    # Set up beam arrays
    models = ('apertif', 'parkes', 'chime', 'gaussian', 'airy', 'apertif_real')
    if model in models:
        place = paths.models() + f'/beams/{model}.npy'
        beam_array = np.load(place)

    # Set up details if using beam arrays
    if model.startswith('apertif'):
        pixel_scale = 0.94/60  # Degrees per pixel [deg]
        beam_size = 25.  # [sq deg]
    elif model == 'parkes':
        pixel_scale = 54/3600  # Degrees per pixel [deg]
        beam_size = 9.  # [sq deg]
    elif model == 'chime':
        pixel_scale = 0.08  # Degrees per pixel [deg]
        beam_size = 180*80  # [sq deg]
        # If integrating over 180*80 degrees
        # Wolfram Alpha input:
        # integral_0^pi integral_0^(4*pi/9) sin(theta) d theta d phi square radians to square degrees
        beam_size = 8522  # [sq deg]
    elif model == 'gaussian':
        pixel_scale = fwhm / 95  # Degrees per pixel [deg]
        beam_size = (pixel_scale*beam_array.shape[0])**2
    elif model == 'airy':
        pixel_scale = fwhm / 31  # Degrees per pixel [deg]
        beam_size = (pixel_scale*beam_array.shape[0])**2
    elif model.startswith('perfect'):
        pixel_scale = None
        beam_size = None
        beam_array = None
    else:
        raise ValueError('Beam model input not recognised.')

    return beam_size, pixel_scale, beam_array


def calc_max_offset(n, fwhm):
    """Calculate the maximum offset of an FRB in an Airy disk.

    Args:
        n (int): Maximum number of wanted sidelobes
        fwhm (float): Size of the FWHM (note it's a diameter) [deg]

    Returns:
        max_offset (float): Maximum offset from beam centre [deg]

    """
    # Allow for cut at FWHM
    if n == 0.5:
        return fwhm / 2  # Max offset = radius of beam

    # Null points of kasin for allow a number of sidelobes
    kasin_nulls = [3.831706, 7.015587, 10.173468, 13.323692, 16.47063,
                   19.615859, 22.760084, 25.903672, 29.046829, 32.18968,
                   35.332308, 38.474766]

    arcsin = np.arcsin(np.deg2rad(fwhm)*kasin_nulls[n]/np.pi)
    if np.isnan(arcsin):
        m = f'Beamsize including {n} sidelobes would be larger than sky \n'
        A = (90/kasin_nulls[n])**2*np.pi
        m += f'Ensure beamsize_at_fwhm is smaller than {A}'
        raise ValueError(m)

    return np.rad2deg(arcsin)


def int_pro_random(shape=(1, 1), pattern='perfect', fwhm=2, max_offset=None,
                   central_freq=1400, beam_array=None, pixel_scale=None):
    """Calculate the intensity profile in random spots of a beam.

    Args:
        shape (tuple): The shape of the array with intensities that you need.

    Returns:
        array, array: intensity profile, offset from beam [degree]

    """
    offset = fwhm/2  # Radius [deg]

    # Take a random location in the 2D beampattern
    offset *= np.sqrt(np.random.random(shape).astype(np.float32))

    # Convert max offset to units of the radius
    if max_offset is not None:
        max_offset /= fwhm/2

    # Allow for a perfect beam pattern in which all is detected
    if pattern.startswith('perfect'):
        offset *= max_offset
        int_pro = np.ones(shape)
        return int_pro, offset

    # Formula's based on 'Interferometry and Synthesis in Radio
    # Astronomy' by A. Richard Thompson, James. M. Moran and
    # George W. Swenson, JR. (Second edition), around p. 15

    if pattern == 'gaussian':
        offset *= max_offset
        alpha = 2*np.sqrt(np.log(2))
        int_pro = np.exp(-(alpha*offset/fwhm)**2)
        return int_pro, offset

    elif pattern == 'airy':
        # Set the maximum offset equal to the null after a sidelobe
        offset *= max_offset
        c = 299792458
        conv = np.pi/180  # Conversion degrees -> radians
        eff_diam = c/(central_freq*1e6*conv*fwhm)
        a = eff_diam/2  # Effective radius of telescope
        lamda = c/(central_freq*1e6)
        ka = (2*np.pi*a/lamda)
        kasin = ka*np.sin(offset*conv)
        int_pro = 4*(j1(kasin)/kasin)**2
        return int_pro, offset

    # Use an array of the beam pattern
    elif beam_array is not None:
        b_shape = beam_array.shape
        ran_x = np.random.randint(0, b_shape[0], shape)
        ran_y = np.random.randint(0, b_shape[1], shape)
        int_pro = beam_array[ran_x, ran_y]
        x_offset = (ran_x-(b_shape[0]/2)) * pixel_scale
        y_offset = (ran_y-(b_shape[1]/2)) * pixel_scale
        offset = go.separation(0, 0, x_offset, y_offset)
        return int_pro, offset

    else:
        raise ValueError(f'Beam pattern "{pattern}" not recognised')


def int_pro_fixed(ra, dec, ra_p, dec_p, lst, pattern='perfect',
                  latitude=0, beam_array=None, pixel_scale=1,
                  mount_type='equatorial'):
    """Calculate intensity profile for fixed location in beam.

    Args:
        ra (array): Right ascension of objects [deg]
        dec (array): Declination of objects [deg]
        ra_p (float): Right ascension of pointing [deg]
        dec_p (float): Declination of pointing [deg]
        lst (float): Local Sidereal Time [deg]

    Returns:
        type: Description of returned object.

    """
    # Weed out perfect beam
    if pattern.startswith('perfect'):
        return np.ones_like(ra), None, None

    # Check input
    args = [ra_p, dec_p, lst, latitude]
    for a in args:
        if a is None:
            raise ValueError('Missing required input')

    # Convert input decimal degrees to radians
    ra = np.deg2rad(ra)
    dec = np.deg2rad(dec)
    ra_p, dec_p, lst, lat = np.deg2rad(args)

    if mount_type == 'equatorial':
        # Convert input coordinates to offset in ra and dec
        dx, dy = go.coord_to_offset(ra_p, dec_p, ra, dec)
    elif mount_type == 'azimuthal':
        # Convert input right ascension to hour angle
        ha = lst - ra
        ha_p = lst - ra_p
        if isinstance(ha, np.ndarray):
            ha[ha > np.pi] -= 2*np.pi
            ha[ha < -np.pi] += 2*np.pi
        else:
            if ha > np.pi:
                ha -= 2*np.pi
            if ha < -np.pi:
                ha += 2*np.pi
        if isinstance(ha_p, np.ndarray):
            ha_p[ha_p > np.pi] -= 2*np.pi
            ha_p[ha_p < -np.pi] += 2*np.pi
        else:
            if ha_p > np.pi:
                ha_p -= 2*np.pi
            if ha_p < -np.pi:
                ha_p += 2*np.pi
        # Convert ha, dec to az, alt
        az, alt = go.hadec_to_azalt(ha, dec, lat)
        az_p, alt_p = go.hadec_to_azalt(ha_p, dec_p, lat)
        # Convert to offset
        # Only valid for +/-30 from the centre
        dx, dy = go.coord_to_offset(az_p, alt_p, az, alt)
    elif mount_type == 'transit':
        # A transit telescope always looks straight up
        az_p = np.deg2rad(180)
        alt_p = np.deg2rad(90)

        # Convert input right ascension to hour angle
        ha = lst - ra
        if isinstance(ha, np.ndarray):
            ha[ha > np.pi] -= 2*np.pi
            ha[ha < -np.pi] += 2*np.pi
        else:
            if ha > np.pi:
                ha -= 2*np.pi
            if ha < -np.pi:
                ha += 2*np.pi

        # Convert ha, dec to az, alt
        az, alt = go.hadec_to_azalt(ha, dec, lat)

        def pol2cart(rho, phi):
            x = rho * np.cos(phi)
            y = rho * np.sin(phi)
            return x, y

        # Convert to offset
        dx, dy = pol2cart(np.pi/2-alt, az+np.pi/2)
    else:
        raise ValueError(f'Invalid mount type: {mount_type}')

    # Convert offsets dx, dy to pixel in beam pattern (round)
    dx_px = (np.round(np.rad2deg(dx) / pixel_scale)).astype(int)
    dy_px = (np.round(np.rad2deg(dy) / pixel_scale)).astype(int)
    ny, nx = beam_array.shape
    x = (nx/2 + dx_px).astype(np.int32)
    y = (ny/2 + dy_px).astype(np.int32)

    # Get the value at this pixel (zero if outside beam pattern)
    if not isinstance(x, np.ndarray):
        x = np.array(x, dtype=np.int32)
        y = np.array(y, dtype=np.int32)

    i, j = beam_array.shape
    outside = ((x < 0) | (x > j-1) | (y < 0) | (y > i-1))

    x[outside] = 0  # Nans don't work in int arrays
    y[outside] = 0

    intensity = np.asarray(beam_array[y, x])
    intensity[((x == 0) & (y == 0))] = np.nan


    return intensity, np.rad2deg(dx), np.rad2deg(dy)
