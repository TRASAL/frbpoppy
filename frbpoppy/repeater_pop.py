"""Allow for repeating FRB sources."""

from frbpoppy.cosmic_pop import CosmicPopulation


class RepeaterPopulation(CosmicPopulation):
    """Allow repeating FRBs to be modeled."""

    def __init__(self,
                 frac=0.05,
                 **kw):
        """Allow repeating FRB sources to be modeled.

        Args:
            frac (float): Fraction of population to be repeaters.
            **kw (all): All of the arguments available to CosmicPopulation.
        """
        super(RepeaterPopulation, self).__init__(**kw)
        self.name = 'repeater'
        self.frac = 0.5
        self.n_rep = round(self.frac*self.n_gen)

    def test(self):
        """Test something."""
        pass


if __name__ == '__main__':
    pop = RepeaterPopulation(n_gen=200, generate=False)
    print(pop.__dict__)
