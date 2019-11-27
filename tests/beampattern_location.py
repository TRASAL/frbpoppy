"""Calculate the location in a beampattern."""
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

from frbpoppy import paths
import frbpoppy.galacticops as go
from convenience import plot_aa_style

pointing = (0, 90)  # [RA, Dec]
point = (0, 89)  # [RA, Dec]
telescope = 'apertif'  # 1 pixel = 0.94'
lon = 6.6
lat = 52
lon = 0
lat = 90
n_pointings = 10
t_obs = (24/n_pointings)*60*60
mount = 'equatorial'
mount = 'altaz'
actual_times = False

place = paths.models() + f'/beams/{telescope}.npy'
beam_array = np.load(place)

if telescope == 'apertif':
    deg_to_pix = 1/(0.94/60)  # [deg]

centre_ij = (np.array(beam_array.shape) / 2).astype(int)

# Generate local sidereal times
step = t_obs/(60*60)*(360/23.9344696)
start = np.random.uniform(0, 360, 1)
lsts = (np.arange(start, 360+start, step) + lon) % 360

if actual_times:
    # Pointings are randomly placed between the year 2000 and 2100
    date_min = go.random_date(datetime(2000, 1, 1), datetime(2100, 1, 1))
    date_max = date_min + timedelta(seconds=int(t_obs*n_pointings))
    time_delta = np.timedelta64(t_obs, 's')
    times = np.arange(date_min, date_max, time_delta, dtype='datetime64')
    lsts = go.datetime_to_gmst(times) + lon


def lst_to_proj(lst, lat=lat):
    """
    Convert HA, Dec to parallactic angle.

    This is the SB rotation w.r.t. the RA-Dec frame
    :param ha: hour angle with unit
    :param dec: declination with unit
    :param lat: Latitude with unit (default: WSRT)
    """
    lat = np.deg2rad(lat)
    ha = np.deg2rad(lst - pointing[0])
    dec = np.deg2rad(pointing[1])
    q = np.arctan2(np.cos(lat)*np.sin(ha),
                  (np.sin(lat)*np.cos(dec)-np.cos(lat)*np.sin(dec)*np.cos(ha)))
    return np.rad2deg(q)


def sph2cart(ra, dec):
    """Spherical coordinates to Cartesian."""
    ra = np.deg2rad(ra)
    dec = np.deg2rad(dec)
    cos_dec = np.cos(dec)
    return np.array([cos_dec*np.cos(ra), cos_dec*np.sin(ra), np.sin(dec)])


def rot(coords, pointing, angle):
    """Rotate coords around pointing with angle.

    Adaptation of Rodrigues' rotation formula
    https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    """
    cs = sph2cart(coords[0], coords[1])
    k = sph2cart(pointing[0], pointing[1])
    q = np.deg2rad(angle)

    sin = np.sin
    cos = np.cos

    rot_cs = cos(q)*cs + sin(q)*np.cross(k, cs)+np.dot(k, cs)*(1-cos(q))*k
    rot_x = rot_cs[0]
    rot_y = rot_cs[1]
    rot_z = rot_cs[2]

    rot_ra = np.arctan2(rot_y, rot_x)
    rot_dec = np.arcsin(rot_z)

    rot_ra = np.rad2deg(rot_ra)
    rot_dec = np.rad2deg(rot_dec)

    return rot_ra, rot_dec


def int_from_pattern(centre, point, lst, test=False):
    """Return intensity value given a pointing and coordinates of an FRB."""
    ra_c, dec_c = centre  # [np.deg2rad(c) for c in centre]
    ra_p, dec_p = point  # [np.deg2rad(p) for p in point]

    # ha = lst - ra_p
    # a = np.sin(dec_p)*np.sin(lat)+np.cos(dec_p)*np.cos(lat)*np.cos(ha)
    # alt = np.arcsin(a)
    # cos_a = (np.sin(dec_p) - np.sin(alt)*np.sin(lat))/(np.cos(alt)*np.cos(lat))
    # az = np.arccos(cos_a)
    # if isinstance(az, np.ndarray):
    #     m = (np.sin(ha) >= 0)
    #     az[m] = 2*np.pi - az[m]
    # else:
    #     if np.sin(ha) >= 0:
    #         az = 2*np.pi - az
    # angle = za

    dec_diff = go.separation(ra_c, dec_c, ra_c, dec_p)
    ra_diff = go.separation(ra_p, dec_p, ra_c, dec_p)
    point_i = centre_ij[1] + int(ra_diff*deg_to_pix)
    point_j = centre_ij[0] + int(dec_diff*deg_to_pix)
    # distance_dec = (dec_p - dec_c)
    # distance_ra = (ra_p - ra_c)*np.cos(dec_p)
    # point_i = centre_ij[1] + int(distance_ra*deg_to_pix)
    # point_j = centre_ij[0] + int(distance_dec*deg_to_pix)

    try:
        int_pro = beam_array[point_i, point_j]

        if test is False:
            return int_pro
        else:
            return int_pro, point_i, point_j
    except IndexError:
        pass


ims = []
for i, lst in enumerate(lsts):

    new_point = False
    q = 0.
    if mount == 'altaz':
        q = lst_to_proj(lst, lat)
        new_point = rot(point, pointing, -q)

    if new_point:
        p = new_point
    else:
        p = point
    a = int_from_pattern(pointing, p, lst, test=True)
    print(f'{point} around {pointing} with angle {q:.3} -> {new_point}, {a}')
    if a:

        beam_array[a[-1], a[-2]] = 1


plot_aa_style()
extent = [-beam_array.shape[0]/2*1/deg_to_pix,
          beam_array.shape[0]/2*1/deg_to_pix,
          -beam_array.shape[1]/2*1/deg_to_pix,
          beam_array.shape[1]/2*1/deg_to_pix]
plt.imshow(beam_array, origin='lower', extent=extent)

plt.savefig(f'./tests/plots/bp.pdf')
