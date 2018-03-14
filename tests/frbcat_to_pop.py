"""Convert a Pandas DataFrame with Frbcat to a Population class."""
from frbpoppy.frbcat import Frbcat
from frbpoppy.do_survey import observe

frbcat = Frbcat().to_pop()
surv_pop = observe(frbcat, 'APERTIF')

frbs = []
fluences = []

for src in surv_pop.sources:
    frbs.append(src.name)
    for frb in src.frbs:
        fluences.append((src.name, frb.snr))

for f in fluences:
    print(f)
