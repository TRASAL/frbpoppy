"""Simulate large CHIME population."""
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from tqdm import tqdm

from frbpoppy import LargePopulation, CosmicPopulation, log10normal, Survey
from frbpoppy import plot, int_pro_fixed, SurveyPopulation

from convenience import plot_aa_style, rel_path


N_SRCS = 1000
N_DAYS = 1
RATE = 1000  # per day

r = CosmicPopulation(N_SRCS, n_days=N_DAYS, repeaters=True)
r.set_dist(model='vol_co', z_max=2.5)
r.set_dm_host(model='gauss', mean=100, std=0)
r.set_dm_igm(model='ioka', slope=1000, std=0)
r.set_dm(mw=False, igm=True, host=True)
r.set_emission_range(low=100e6, high=10e9)
r.set_lum(model='powerlaw', per_source='different', low=1e40, high=1e45,
          power=0)
r.set_si(model='gauss', mean=0, std=0)
r.set_w(model='log10normal', per_source='different', mean=0.1, std=1)
# rate = log10normal(RATE, 1, int(N_SRCS))
rate = RATE
r.set_time(model='poisson', rate=rate)

# Set up survey
s = Survey('chime', n_days=N_DAYS)
s.set_beam(model='chime')


# Only generate FRBs in CHIME's survey region
r.set_direction(model='uniform',
                min_ra=s.ra_min,
                max_ra=s.ra_max,
                min_dec=s.dec_min,
                max_dec=s.dec_max)

r.generate()

surv_pop = SurveyPopulation(r, s, test_beam_placement=True)


#
# surv_pop = LargePopulation(r, s, run=True).pops[0]
# surv_pop.name = 'chime'
# surv_pop.save()

# # import IPython; IPython.embed()
print(surv_pop.source_rate)
print(surv_pop.burst_rate)

# Set up image properties
plot_aa_style()
fig, ax = plt.subplots()

# Get beam properties
beam_array = s.beam_array
pattern = s.beam_pattern
pixel_scale = np.float64(s.pixel_scale)
latitude = s.latitude
mount_type = s.mount_type

# Plot image in log scale
# Set up beam pattern array to plot
ny, nx = beam_array.shape
xx = np.arange(-nx/2, nx/2) * pixel_scale
yy = np.arange(-ny/2, ny/2) * pixel_scale

# Find corresponding Ra, Dec coordinates
extent = np.array([xx[0], xx[-1], yy[0], yy[-1]])
norm = LogNorm(vmin=np.min(beam_array), vmax=np.max(beam_array))
img = plt.imshow(beam_array, extent=extent, norm=norm, aspect='auto')

# Set up a tuple of pointings if not given
dx, dy = surv_pop.dxys

ax.scatter(np.array(dx), np.array(dy), color='w', marker='o', s=(72/300)**2)

# Set axis labels
ax.set_xlabel(r'X Offset ($^{\circ}$)')
ax.set_ylabel(r'Y Offset ($^{\circ}$)')

# Set axis limits
ax.set_xlim(extent[0], extent[1])
ax.set_ylim(extent[2], extent[3])

# Add colour bar
cb = plt.colorbar(img)
cb.set_label('Intensity')

# Save figure
plt.tight_layout()
plt.savefig(rel_path('./plots/det_in_chime_beam.pdf'))
