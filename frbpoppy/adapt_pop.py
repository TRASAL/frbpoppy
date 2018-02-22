"""Set of functions which adapt a population."""

import random

import frbpoppy.galacticops as go
import frbpoppy.distributions as dis


class Adapt:
    """Adapt a population."""

    def __init__(self, pop):
        """Set the population."""
        self.pop = pop

    def dm_host(self, dm):
        """
        Adapt the default dispersion measure of an frb.

        Args:
            dm (int): Dispersion measure host

        Returns:
            pop (Population): Population with new parameters

        """
        for src in self.pop.sources:
            src.dm_host = dm
            src.dm = src.dm_mw + src.dm_igm + src.dm_host

        # Save the population
        self.pop.pickle_pop()
        return self.pop

    def dm_igm(self, slope):
        """
        Adapt the slope of the intergalactic dispersion measure.

        Args:
            slope (int): Dispersion measure slope

        Returns:
            pop (Population): Population with new parameters

        """
        for src in self.pop.sources:
            src.dm_igm = go.ioka_dm_igm(src.z, slope=slope)
            src.dm = src.dm_mw + src.dm_igm + src.dm_host

        # Save the population
        self.pop.pickle_pop()
        return self.pop

    def freq(self, freq_min, freq_max):
        """
        Adapt the range of frequencies over which an frb emits.

        Args:
            freq_min (int): Minimum frequency
            freq_max (int): Maximum frequency

        Returns:
            pop (Population): Population with new parameters

        """
        self.pop.f_min = freq_min
        self.pop.f_max = freq_max
        self.pop.pickle_pop()
        return self.pop

    def lum_bol(self, lum_min, lum_max, lum_pow):
        """
        Adapt the range of luminosities over which an frb can be emitted.

        Args:
            lum_min (int): Minimum luminosity
            lum_max (int): Maximum luminosity
            lum_pow (int): Power of the luminosity slope

        Returns:
            pop (Population): Population with new parameters

        """
        for source in self.pop.sources:
            for frb in source.frbs:
                frb.lum_bol = dis.powerlaw(lum_min, lum_max, lum_pow)

        # Save the population
        self.pop.pickle_pop()
        return self.pop

    def si(self, si_mean, si_sigma):
        """
        Adapt the population's mean spectral index and the standard deviation.

        Args:
            si_mean (float): Mean spectral index
            si_sigma (float): Mean standard deviation of the spectral index

        Returns:
            pop (Population): Population with new parameters

        """
        for source in self.pop.sources:
            for frb in source.frbs:
                frb.si = random.gauss(self.pop.si_mean, self.pop.si_sigma)

        # Save the population
        self.pop.pickle_pop()
        return self.pop

    def w_int(self, w_int_min, w_int_max):
        """
        Adapt the pulse width of a population.

        Args:
            w_int_min (float): Minimum pulse width [ms]
            w_int_max (float): Maximum pulse width [ms]

        Returns:
            pop (Population): Population with new parameters

        """
        # Adapt just the intrinsic pulse width
        for source in self.pop.sources:
            for frb in source.frbs:
                frb.w_int = dis.redshift_w(z=source.z,
                                           w_min=w_int_min,
                                           w_max=w_int_max)

        # Save the population
        self.pop.pickle_pop()
        return self.pop
