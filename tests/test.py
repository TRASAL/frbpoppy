"""Test calculating transit times for location on Earth."""
import numpy as np
import matplotlib.pyplot as plt

latitude_obs = 0.  # deg


def calc_transit_time(lat, dec, unit='deg', beamsize=20626.5):
    """Calculate total time of an object above horizon.

    Args:
        lat (float): Latitude of observatory in radians.
        dec (array): Declination of celestial object in radians.
        unit (str): 'deg' or 'rad'
        beamsize (float): Beam size in sq. deg

    Returns:
        float: Fraction of time above horizon in fractional days

    """
    if unit == 'deg':
        lat = np.deg2rad(lat)
        dec = np.deg2rad(dec)

    # Beamsize to radius in fractional degrees
    r = np.sqrt(beamsize/np.pi)
    extra_lim = np.deg2rad(90 - r)

    times = np.ones_like(dec)
    lim = np.pi/2. - lat
    always_visible = dec > lim
    never_visible = dec < -lim
    sometimes_visible = ((dec > -lim) & (dec < lim))

    sm = sometimes_visible
    times[sm] = 2*np.rad2deg(np.arccos(-np.tan(dec[sm])*np.tan(lat)))/360.

    times[always_visible] = 1.
    times[never_visible] = 0.

    return times


decs = np.linspace(-90, 90, 100)
lat = latitude_obs
times = calc_transit_time(lat, decs)
plt.ylabel('Fraction of day at which visible')
plt.xlabel('Declination sources')
plt.plot(decs, times)
plt.show()
