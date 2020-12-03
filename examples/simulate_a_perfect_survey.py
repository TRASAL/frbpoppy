"""Example of simulating a perfect survey."""
from frbpoppy import CosmicPopulation, Survey, SurveyPopulation, plot

# Generate an FRB population
cosmic_pop = CosmicPopulation.simple(1e4, generate=True)

# Setup a survey
survey = Survey('perfect')

# Observe the FRB population
survey_pop = SurveyPopulation(cosmic_pop, survey)

# Check the detection rates
print(survey_pop.source_rate)

# Note that due to redshift you won't see all bursts, as some will have
# redshifted out of the observing time. But no matter how faint, you'll see
# all bursts within the observing time

# Plot populations
plot(cosmic_pop, survey_pop, tns=False, mute=False)
