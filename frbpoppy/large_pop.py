"""Methods to deal with large populations."""
import numpy as np
import os
import uuid
from tqdm import tqdm

from frbpoppy.misc import pprint
from frbpoppy.survey_pop import SurveyPopulation
from frbpoppy.population import unpickle
from frbpoppy.paths import paths


class LargePopulation:
    """Method to deal with a large population.

    Basically splits it up, and merges the results.
    """

    def __init__(self, pop, *surveys, max_size=1e6, run=True):
        """Set arguments."""
        self.pop = pop
        self.base_name = pop.name
        self.surveys = surveys
        self.max_size = int(max_size)
        self.uids = None

        if run:
            self.run()
            self.merge()

    def run(self):
        """Run the generating and surveying of a large population."""
        pprint(f'Running a large {self.base_name} population')
        d = divmod(self.pop.n_srcs, self.max_size)
        sizes = [self.max_size for i in range(d[0])]
        if d[1] != 0:
            sizes.append(d[1])
        self.uids = [str(uuid.uuid4())[:8] for s in sizes]

        for i, n in enumerate(tqdm(sizes, desc='Subpopulations')):
            pop = self.pop
            pop.n_srcs = n
            pop.uid = self.uids[i]
            pop.generate()

            for surv in self.surveys:
                surv_pop = SurveyPopulation(pop, surv)
                surv_pop.name = f'{self.base_name}_{surv_pop.name}'
                surv_pop.uid = pop.uid
                surv_pop.save()

    def merge(self):
        """Merge populations."""
        pprint('Merging populations')

        self.pops = []

        for s in self.surveys:
            pops = []
            files = [f'{self.base_name}_{s.name}_{uid}' for uid in self.uids]
            for f in files:
                pops.append(unpickle(f))

            # Main population
            mp = pops[0]

            # Merge each parameter
            for attr in mp.frbs.__dict__.keys():
                parm = getattr(mp.frbs, attr)
                if type(parm) is np.ndarray:
                    parms = []
                    for pop in pops:
                        parms.append(getattr(pop.frbs, attr))

                    try:
                        merged_parm = np.concatenate(parms, axis=0)
                    except ValueError as e:
                        # Check maximum size values should be padded to
                        max_size = max([p.shape[1] for p in parms])
                        new_parms = []

                        # Ensure matrices are the same shapes by padding them
                        for p in parms:
                            if p.shape[1] != max_size:
                                padded_p = np.zeros((p.shape[0], max_size))
                                padded_p[:] = np.nan
                                padded_p[:, :p.shape[1]] = p
                                new_parms.append(padded_p)
                            else:
                                new_parms.append(p)

                        merged_parm = np.concatenate(new_parms, axis=0)

                    setattr(mp.frbs, attr, merged_parm)

            # Add up detections
            mp.rate.det = len(mp.frbs.snr)
            for pop in pops[1:]:
                mp.rate.faint += pop.rate.faint
                mp.rate.late += pop.rate.late
                mp.rate.out += pop.rate.out

            # Recalculate detection rates
            mp.calc_rates(s)

            # Save the main population as one big population
            mp.uid = None
            mp.save()

            # Remove all of the smaller ones
            for f in files:
                p = paths.populations() + f'{f}.p'
                os.remove(p)

            self.pops.append(mp)


def main():
    """Compare a standard Survey Population with a Large Survey Population."""
    from frbpoppy.cosmic_pop import CosmicPopulation
    from frbpoppy.survey import Survey
    from frbpoppy.survey_pop import SurveyPopulation
    from frbpoppy.do_plot import plot

    pop = CosmicPopulation(int(5e5), name='test', generate=False)

    surveys = []
    for s in ['apertif']:
        surveys.append(Survey(s))

    large_pops = LargePopulation(pop, *surveys).pops

    # For comparison
    surv_pop = SurveyPopulation(pop, surveys[0])

    plot(*large_pops, surv_pop, mute=False, frbcat=False)


if __name__ == '__main__':
    main()
