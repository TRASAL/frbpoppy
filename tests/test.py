@profile
def main():

    import numpy as np
    from frbpoppy import CosmicPopulation, Survey, SurveyPopulation

    pop = CosmicPopulation(int(1e4))
    s = Survey('apertif')
    surv_pop = SurveyPopulation(pop, s)
    del pop
    tot = 0

    for attr in surv_pop.frbs.__dict__.keys():
        parm = getattr(surv_pop.frbs, attr)
        if type(parm) is np.ndarray:
            print(attr, parm.dtype, parm.nbytes)
            tot += parm.nbytes

    print(tot)


if __name__ == '__main__':
    main()
