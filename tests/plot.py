from frbpoppy import plot

surveys = ('htru', 'fast', 'puma-full', 'chord', 'ska1-low', 'ska1-mid')
pops = [f'future_complex_{s}' for s in surveys]
plot(*pops, frbcat=False, mute=False)
