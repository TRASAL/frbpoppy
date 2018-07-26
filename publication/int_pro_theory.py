"""Plot intensity profile of theoretical beam patterns."""
import matplotlib.pyplot as plt
from frbpoppy.survey import Survey

PATTERNS = ['perfect', 'tophat', 'gaussian', 'airy']
SURVEY = 'HTRU'
MIN_Y = 1e-7
n = 50000

for pattern in PATTERNS:

    s = Survey(SURVEY, gain_pattern=pattern)

    if pattern != 'tophat':
        data = [s.intensity_profile(test=True) for e in range(n)]

        # Clean up lower limit
        data = [d for d in data if d[1] >= MIN_Y]
        offset, int_pro = zip(*sorted(data))

        # Offset in degrees
        offset = [o/60. for o in offset]

        plt.plot(offset, int_pro, label=pattern)

    # Make a nicer looking tophat
    else:
        x_min, x_max = min(offset), max(offset)
        midway = (x_max - x_min)/2 + x_min
        plt.plot([x_min, midway, midway], [1, 1, MIN_Y], label='tophat')


plt.xlabel(f'Offset ($\degree$)')
plt.ylabel('Intensity Profile')
plt.yscale('log')
plt.legend()
plt.tight_layout()
# plt.show()
plt.savefig('plots/int_pro_theory.pdf')
