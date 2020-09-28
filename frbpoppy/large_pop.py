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
                surv_pop = SurveyPopulation(pop, surv, scale_by_area=False)
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
                    except ValueError:
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
            for pop in pops[1:]:
                mp.source_rate.faint += pop.source_rate.faint
                mp.source_rate.late += pop.source_rate.late
                mp.source_rate.out += pop.source_rate.out
                mp.source_rate.det += pop.source_rate.det

                if mp.repeaters:
                    mp.burst_rate.faint += pop.burst_rate.faint
                    mp.burst_rate.late += pop.burst_rate.late
                    mp.burst_rate.out += pop.burst_rate.out
                    mp.burst_rate.det += pop.burst_rate.det
                    mp.burst_rate.pointing += pop.burst_rate.pointing

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

    pop = CosmicPopulation.simple(int(3e3), generate=False)
    pop.name = 'large'
    survey = Survey('perfect')
    large_pop = LargePopulation(pop, survey, max_size=1e3).pops[0]

    # For comparison
    surv_pop = SurveyPopulation(pop, survey)
    surv_pop.name = 'normal'

    plot(large_pop, surv_pop, mute=False, tns=False)


if __name__ == '__main__':
    main()
