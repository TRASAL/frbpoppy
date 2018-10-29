import matplotlib.pyplot as plt
from astropy.cosmology import Planck15 as cosmo
import numpy as np

z = np.logspace(-10, 1, 10000)
dc = cosmo.comoving_distance(z).value
f = 1/((1+z)*dc**2)
d = z.copy()

vol = cosmo.comoving_volume(z)
vol_bins = np.random.uniform(min(vol.value), max(vol.value), 10000)
z_bin = z[np.where(vol_bins[-1] <= vol.value)[0][0]]


def co_vol_to_z(b):
    return z[np.where(b <= vol.value)[0][0]]


co_vol_to_z = np.vectorize(co_vol_to_z)
z_bins = co_vol_to_z(vol_bins)
dc_bins = cosmo.comoving_distance(z_bins).value
f_bins = 1/((1+z_bins)*dc_bins**2)
d_bins = z_bins.copy()

# plt.hist(f, histtype='step', log=True, label='Original')
# plt.hist(f_bins, histtype='step', log=True, label='Comoving Volume')

plt.plot(z, f, '.')
plt.yscale('log')
plt.xlabel('z')
plt.ylabel('Fluence')
# plt.xscale('log')
plt.legend()
