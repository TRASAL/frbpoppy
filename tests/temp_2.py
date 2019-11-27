"""Rotate coordinates."""
import numpy as np


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


coord = (90, 0)
angle = 85
pointing = (0, 0)

ra, dec = rot(coord, pointing, angle)
print(f'Rotate {coord} around {pointing} by {angle} -> {ra}, {dec}')
