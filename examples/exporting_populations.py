"""Examples of saving populations."""
from frbpoppy import CosmicPopulation, unpickle

# Set up an FRB population
pop = CosmicPopulation.simple(1e4, generate=True)

# Set filename
file_name = 'saving_example'

# Exporting a population as a pickled object
pop.to_pickle(f'./{file_name}.p')

# Unpickling a population
copy_of_pop = unpickle(f'./{file_name}.p')

# Exporting frb information to csv
pop.frbs.to_csv(f'./{file_name}.csv')

# Exporting frb information to a Pandas DataFrame
df = pop.frbs.to_df()

# Alternatively save in '/data/results' as a pickled file
pop.name = file_name
pop.save()
