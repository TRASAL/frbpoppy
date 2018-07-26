"""Plot intensity profile of sidelobes."""
import matplotlib.pyplot as plt
from frbpoppy.survey import Survey

SIDELOBES = [0, 1, 2, 8]
SURVEY = 'APERTIF'
MIN_Y = 1e-7
n = 50000

for sidelobe in reversed(SIDELOBES):

    args = {'sidelobes': sidelobe,
            'equal_area': True}

    s = Survey(SURVEY, gain_pattern='airy', **args)

    data = [s.intensity_profile(test=True, **args) for e in range(n)]

    # Clean up lower limit
    data = [d for d in data if d[1] >= MIN_Y]
    offset, int_pro = zip(*sorted(data))

    # Offset in degrees
    offset = [o/60. for o in offset]

    label = f'{sidelobe} sidelobes'
    if sidelobe == 1:
        label = label[:-1]

    plt.plot(offset, int_pro, label=label)

plt.xlabel(f'Offset ($\degree$)')
plt.ylabel('Intensity Profile')
plt.yscale('log')
plt.legend()
plt.tight_layout()
# plt.show()
plt.savefig('plots/int_pro_sidelobes.pdf')
