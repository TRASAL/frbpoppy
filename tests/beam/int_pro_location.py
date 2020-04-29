"""Plot frb locations in an intensity profile."""
import matplotlib.pyplot as plt
import numpy as np

from frbpoppy import Survey
from convenience import plot_aa_style, rel_path

surv = 'parkes'

survey = Survey(surv)
survey.beam_pattern = 'gaussian'

# Set up some testing parameters
ra = np.array([181, 182])
dec = np.array([89, 88.5])
ra_pointing = 181.
dec_pointing = 89
lst = 90.

# Set up beam properties
_, _ = survey.calc_beam(shape=0)

# Calculate position in beam
int, x, y = survey.calc_fixed_int_pro(ra=ra,
                                      dec=dec,
                                      ra_p=ra_pointing,
                                      dec_p=dec_pointing,
                                      lst=lst,
                                      test=True)


plot_aa_style()
fig, ax = plt.subplots()

# Plot beam
ny, nx = survey.beam_array.shape
xx = np.arange(-nx/2, nx/2) * survey.pixel_scale   # [deg]
yy = np.arange(-ny/2, ny/2) * survey.pixel_scale   # [deg]
# Rasterise to reduce image size
img = plt.pcolormesh(*np.meshgrid(xx, yy), survey.beam_array, rasterized=True)

# Plot position in beam
try:
    plt.scatter(xx[x], yy[y], marker='x', color='r')
except IndexError:  # Outside of beam
    pass

# Plot properties
plt.xlabel(r'Offset ($^{\circ}$)')
plt.ylabel(r'Offset ($^{\circ}$)')
ax.set_aspect('equal')
plt.tight_layout()
plt.savefig(rel_path('./plots/int_pro_location.pdf'), dpi=300)
