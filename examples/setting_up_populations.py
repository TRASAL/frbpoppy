"""Setting up various cosmic populations."""
from frbpoppy import CosmicPopulation

# You can set up a population with arguments ...
pop = CosmicPopulation(1e4, n_days=1, name='my_own_population', repeaters=True,
                       generate=False)
# ... but also adapt specific components:

# The numer density / distance parameters
pop.set_dist(model='vol_co', z_max=0.01, alpha=-1.5,
             H_0=67.74, W_m=0.3089, W_v=0.6911)

# Which dispersion measure components to include
pop.set_dm(mw=True, igm=True, host=True)

# Dispersion measure properties
pop.set_dm_host(model='gauss', mean=100, std=200)
pop.set_dm_igm(model='ioka', slope=1000, std=None)
pop.set_dm_mw(model='ne2001')

# Emission range of FRB sources
pop.set_emission_range(low=100e6, high=10e9)

# Luminsity of FRBs
# See the per_source argument? That always you to give different properties
# to different bursts from the same source. You can do that for the luminosity,
# or any of the following parameters
pop.set_lum(model='powerlaw', low=1e38, high=1e38, power=0,
            per_source='different')

# Pulse width
pop.set_w(model='uniform', low=10, high=10)

# Spectral index
pop.set_si(model='gauss', mean=0, std=0)

# If repeaters, how they repeat
pop.set_time(model='regular', rate=2)

# And then generate the population!
pop.generate()

# Or simply use some predefined models
pop_simple = CosmicPopulation.simple(1e4, generate=True)
pop_complex = CosmicPopulation.complex(1e4, generate=True)
# Note that for generating these you have to set generate=True
