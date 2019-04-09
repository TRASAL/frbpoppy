"""Use standard populations to speed up calculation times."""
import os
from frbpoppy import CosmicPopulation, SurveyPopulation, paths, unpickle


class StandardCosmicPops:
    """docstring for StandardCosmicPop."""

    def __init__(self, sort, size, alpha, gamma):
        """Quickly get standard populations.

        Args:
            sort (str): Which type of population, standard, std_candle etc.
            size (str): Choice of 'small', 'medium' or 'large'

        Returns:
            pop: Desired population.

        """
        self.sort = sort
        self.size = size
        self.alpha = alpha
        self.gamma = gamma
        self.name = f'{self.sort}_{self.size}'

        if alpha is not None:
            self.name += f'_{self.alpha}'

        if gamma is not None:
            self.name += f'_{self.gamma}'

        self.path = paths.populations() + self.name + '.p'

        if self.size == 'small':  # For testing purposes
            self.n = int(1e5)
        elif self.size == 'medium':  # For laptop purposes
            self.n = int(1e6)
        elif self.size == 'large':  # For cluster purposes
            self.n = int(1e7)
        elif self.size == 'huge':  # For cluster purposes
            self.n = int(1e8)

    def standard_pop(self):
        """Generate a standard population."""
        pop = CosmicPopulation(self.n,
                               days=1,
                               name=self.name,
                               H_0=67.74,
                               W_m=0.3089,
                               W_v=0.6911,
                               dm_host_model='normal',
                               dm_host_mu=100,
                               dm_host_sigma=0,
                               dm_igm_index=1000,
                               dm_igm_sigma=None,
                               dm_mw_model='ne2001',
                               emission_range=[10e6, 10e9],
                               lum_range=[1e40, 1e45],
                               lum_index=0.,
                               n_model='vol_co',
                               alpha=-1.5,
                               pulse_model='uniform',
                               pulse_range=[1., 1.],
                               pulse_mu=1.6,
                               pulse_sigma=1.,
                               si_mu=-1.4,
                               si_sigma=1.,
                               z_max=2.5)

        pop.save()
        return pop

    def standard_candle_pop(self):
        """Generate a standard candle population."""
        pop = CosmicPopulation(self.n,
                               days=1,
                               name=self.name,
                               H_0=67.74,
                               W_m=0.3089,
                               W_v=0.6911,
                               dm_host_model='normal',
                               dm_host_mu=100,
                               dm_host_sigma=0,
                               dm_igm_index=1000,
                               dm_igm_sigma=None,
                               dm_mw_model='ne2001',
                               emission_range=[10e6, 10e9],
                               lum_range=[1e36, 1e36],
                               lum_index=0.,
                               n_model='sfr',
                               alpha=-1.5,
                               pulse_model='uniform',
                               pulse_range=[1., 1.],
                               pulse_mu=1.6,
                               pulse_sigma=0.,
                               si_mu=0.,
                               si_sigma=0.,
                               z_max=2.5)

        pop.save()
        return pop

    def alpha_pop(self):
        """Generate a population varying with alpha."""
        pop = CosmicPopulation(self.n,
                               days=1,
                               name=self.name,
                               H_0=67.74,
                               W_m=0.3089,
                               W_v=0.6911,
                               dm_host_model='normal',
                               dm_host_mu=100,
                               dm_host_sigma=0,
                               dm_igm_index=1000,
                               dm_igm_sigma=None,
                               dm_mw_model='ne2001',
                               emission_range=[10e6, 10e9],
                               lum_range=[1e40, 1e45],
                               lum_index=0.,
                               n_model='vol_co',
                               alpha=self.alpha,
                               pulse_model='uniform',
                               pulse_range=[1., 1.],
                               pulse_mu=1.6,
                               pulse_sigma=1.,
                               si_mu=-1.4,
                               si_sigma=1.,
                               z_max=2.5)

        pop.save()
        return pop

    def alpha_simple_pop(self):
        """Generate a simple local population varying with alpha."""
        pop = CosmicPopulation(self.n,
                               days=1,
                               name=self.name,
                               H_0=67.74,
                               W_m=0.3089,
                               W_v=0.6911,
                               dm_host_model='normal',
                               dm_host_mu=0.,
                               dm_host_sigma=0.,
                               dm_igm_index=0.,
                               dm_igm_sigma=None,
                               dm_mw_model='zero',
                               emission_range=[10e6, 10e9],
                               lum_range=[1e39, 1e39],
                               lum_index=0.,
                               n_model='vol_co',
                               alpha=self.alpha,
                               pulse_model='uniform',
                               pulse_range=[0.01, 0.01],
                               pulse_mu=1.6,
                               pulse_sigma=1.,
                               si_mu=0.,
                               si_sigma=0.,
                               z_max=0.01)

        pop.save()
        return pop

    def gamma_pop(self):
        """Generate a population varying with spectral index."""
        pop = CosmicPopulation(self.n,
                               days=1,
                               name=self.name,
                               H_0=67.74,
                               W_m=0.3089,
                               W_v=0.6911,
                               dm_host_model='normal',
                               dm_host_mu=0.,
                               dm_host_sigma=0.,
                               dm_igm_index=0.,
                               dm_igm_sigma=None,
                               dm_mw_model='zero',
                               emission_range=[10e6, 10e9],
                               lum_range=[10**44.5, 10**44.5],
                               lum_index=0.,
                               n_model='vol_co',
                               alpha=-1.5,
                               pulse_model='uniform',
                               pulse_range=[0.1, 0.1],
                               pulse_mu=1.6,
                               pulse_sigma=1.,
                               si_mu=self.gamma,
                               si_sigma=0.,
                               z_max=2.5)

        pop.save()
        return pop


def get_cosmic_pop(sort, size, load=True, overwrite=False,
                   alpha=None, gamma=None):
    """Quickly get standard populations.

    Args:
        sort (str): Which type of population, standard, std_candle etc.
        size (str): Choice of 'small', 'medium' or 'large'
        load (bool): Whether to load in a population
        overwrite (bool): Check whether a population has already
            been run. If overwrite is true, it will always make a new
            instance.

    Returns:
        pop: Desired population.

    """
    pop = StandardCosmicPops(sort, size, alpha=alpha, gamma=gamma)

    # Skip loading a population if you don't have to
    if not load:
        return pop.name

    # Go for an earlier version if available
    if not overwrite:
        if os.path.isfile(pop.path):
            return unpickle(pop.path)

    # Else generate a standard population
    if pop.sort == 'standard':
        return pop.standard_pop()
    if pop.sort == 'standard_candle':
        return pop.standard_candle_pop()
    if pop.sort == 'alpha':
        return pop.alpha_pop()
    if pop.sort == 'gamma':
        return pop.gamma_pop()
    if pop.sort == 'alpha_simple':
        return pop.alpha_simple_pop()


def get_survey_pop(pop, survey, overwrite=False):
    """Quickly get survey populations.

    Args:
        pop (CosmicPopulation): Population to survey
        survey (Survey): Survey to use
        overwrite (bool): Check whether a population has already
            been run. If overwrite is true, it will always make a new
            instance.

    Returns:
        pop: Desired population.

    """
    observe = True

    # Check where a possible population would be located
    path = ''
    if isinstance(pop, str):
        name = f'{pop}_{survey.name}'
        path = paths.populations() + name + '.p'

    # Check if the file exists
    if not overwrite:
        if os.path.isfile(path):
            observe = False
            return unpickle(path)

    # If all else fails observe again
    if observe:

        if isinstance(pop, str):
            m = 'No survey population exists, yet no surveying requested'
            raise ValueError(m)

        surv_pop = SurveyPopulation(pop, survey)
        surv_pop.name = f'{pop.name}_{survey.name}'
        surv_pop.save()
        return surv_pop
