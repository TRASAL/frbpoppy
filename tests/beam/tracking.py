"""Plot objects through pointing in beam pattern."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D

from frbpoppy import Survey, int_pro_fixed

from tests.convenience import plot_aa_style, rel_path

T_OBS = 60*60  # seconds
OBJ_OFFSET = 1  # deg
INTENSITY_LIMIT = 1e-8

# Set up survey
survey = Survey('chime', n_days=0.99)
survey.set_beam('chime')

# Get beam properties
beam_array = survey.beam_array
pattern = survey.beam_pattern
pixel_scale = np.float64(survey.pixel_scale)
latitude = survey.latitude
mount_type = survey.mount_type


# Limit extent of intensity for plotting reasons
beam_array[beam_array < INTENSITY_LIMIT] = INTENSITY_LIMIT

# Set pointings
survey.gen_pointings()
pointings = survey.pointings

# Targets
ra = 0
decs = [-25, 25, 75]

# Set up beam pattern array to plot
ny, nx = beam_array.shape
xx = np.arange(-nx/2, nx/2) * pixel_scale
yy = np.arange(-ny/2, ny/2) * pixel_scale

# Find corresponding Ra, Dec coordinates
extent = np.array([xx[0], xx[-1], yy[0], yy[-1]])

# Set up image properties
plot_aa_style()
fig, ax = plt.subplots()

# Plot image in log scale
norm = LogNorm(vmin=np.min(beam_array), vmax=np.max(beam_array))
img = plt.imshow(beam_array, extent=extent, norm=norm, aspect='auto')

# Generate local sidereal times
max_t = survey.n_days
t_obs = survey.t_obs/60/60/24  # days
if T_OBS:
    t_obs = T_OBS/60/60/24  # days
times = np.arange(0, max_t+t_obs, t_obs)  # [days]
lsts = times[:-1]*360*(24/23.9344696) % 360  # Local sidereal time [deg]
# lsts += np.random.uniform(0, 360)  # Add random offset
lsts %= 360

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# Show the tracking of objects over time
# Loop over targets
for i, dec in enumerate(decs):
    # Plot position over time
    for ii, lst in enumerate(lsts):
        ra_p = pointings[0][ii]
        dec_p = pointings[1][ii]

        int_pro, dx, dy = int_pro_fixed(ra, dec,
                                        ra_p, dec_p,
                                        lst,
                                        pattern=pattern,
                                        latitude=latitude,
                                        beam_array=beam_array,
                                        pixel_scale=pixel_scale,
                                        mount_type=mount_type)

        ax.scatter(dx, dy, color=colors[i % 9], marker='x', s=20)

# Set axis labels
ax.set_xlabel(r'East-West Offset ($^{\circ}$)')
ax.set_ylabel(r'North-South Offset ($^{\circ}$)')

# Set axis limits
ax.set_xlim(extent[0], extent[1])
ax.set_ylim(extent[2], extent[3])

# Add colour bar
cb = plt.colorbar(img)
cb.set_label('Intensity')

# Add legend
# Add line styles
elements = []
for i, dec in enumerate(decs):
    color = colors[i % 9]
    label = str(int(dec)) + r' $^{\circ}$'
    elements.append((Line2D([0], [0], marker='x', color=color,
                     linestyle='None'), label))
lines, labels = zip(*elements)
ax.legend(lines, labels, loc='upper center', ncol=3, framealpha=1,
          prop={'size': 8}, bbox_to_anchor=(0.5, 1.07),
          bbox_transform=ax.transAxes)

# Save figure
plt.tight_layout()
plt.savefig(rel_path('./plots/tracking.pdf'), dpi=600)
