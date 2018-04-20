"""Calculate where null points lie."""
from collections import deque

from frbpoppy.survey import Survey

surv = Survey('APERTIF')
surv.gain_pattern = 'airy'

m = 0.01
past = deque([1.0, 1.0], maxlen=2)

for x in range(1000):
    int_pro = surv.intensity_profile(mul=m)
    # print(m, int_pro)

    if int_pro > past[1] and past[1] < past[0]:
        null = m - 0.01
        print(f'Null point found using a multiplication factor of {null:.3}')

    # Reset for new round
    past.append(int_pro)
    m += 0.01
