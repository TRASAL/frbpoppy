import healpy as hp
from frbpoppy.paths import paths
import os


def dm_mw(gl, gb, model='ne2001'):
    """ Return DM values for given coordinates from NE2001 healpix

    Args:
        gl (array): Galactic longitude [fractional degrees]
        gb (array): Galactic latitude [fractional degrees]
        model (str): one of 'ne2001', 'ymw16'

    Returns:
        dm_mw (array): Galactic dispersion measure [pc*cm^-3]
    """
    if model == 'ne2001':
        data_file = os.path.join(paths.models(), 'healpix/dm-ne2001-30kpc.fits')
    elif model == 'ymw16':
        data_file = os.path.join(paths.models(), 'healpix/dm-ymw16-30kpc.fits')
    data = hp.read_map(data_file)
    nside = hp.npix2nside(len(data))  # This sets map resolution  (NSIDE)
    pixloc = hp.ang2pix(nside, gl, gb, lonlat=True)
    return data[pixloc]
