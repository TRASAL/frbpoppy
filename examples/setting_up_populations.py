"""Setting up various cosmic populations."""
from frbpoppy import CosmicPopulation

# You can either set up a population with arguments ...
pop_arg = CosmicPopulation(1e4)
# .. such as which components to use for calculating the dispersion measure
pop_arg.set_dm(mw=False, igm=True, host=False)
pop_arg.generate()

# ... or use some predefined models
pop_simple = CosmicPopulation.simple(1e4, generate=True)
pop_complex = CosmicPopulation.complex(1e4, generate=True)
# Note that for generating these you have to set generate=True
