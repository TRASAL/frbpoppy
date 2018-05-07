"""Calculate where null points of an Airy pattern lie."""

USING_FRBPOPPY = False
STEPSIZE = 1e-6

if USING_FRBPOPPY:
    from collections import deque

    from frbpoppy.survey import Survey

    surv = Survey('FAST')
    surv.gain_pattern = 'airy'

    m = STEPSIZE
    past = deque([1.0, 1.0], maxlen=2)

    print(f'Null points found using a multiplication factors of')

    for x in range(int(1/STEPSIZE)**2):
        int_pro = surv.intensity_profile(mul=m)

        if int_pro > past[1] and past[1] < past[0]:
            null = m - STEPSIZE
            print(f'{null},')

        # Reset for new round
        past.append(int_pro)
        m += STEPSIZE

else:
    from scipy.special import j1
    from numpy import arange
    from bokeh.plotting import figure, show

    x_range = arange(0, 100, STEPSIZE)
    y_range = 4*(j1(x_range)/x_range)**2
    nulls = []

    for i in range(2, len(x_range)):
        if y_range[i-2] > y_range[i-1] < y_range[i]:
            null = (x_range[i-1], y_range[i-1])
            nulls.append(null)

    x_null, y_null = zip(*nulls)
    print(x_null)

    # p = figure(title="Bessel function over x")
    # p.line(x_range, y_range)
    # p.circle(x_null, y_null)
    # show(p)
