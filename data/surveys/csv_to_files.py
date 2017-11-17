import pandas as pd
import os

path = os.path.dirname(os.path.abspath(__file__))
df = pd.read_csv(path + '/surveys.csv',
                 index_col='Parameters')

for c in df:

    print('Writing {} survey parameters'.format(c))

    # Header
    text = ['#'*78]
    text.append('# {} survey parameters'.format(c))
    text.append('# WARNING: Values must be taken with a grain of salt - guesses abound!')
    text.append('#'*78)

    # Parameters
    for i, r in df[c].iteritems():
        if r != '#':
            text.append('{:8} ! {}'.format(r, i))

    # Output
    with open(path + '/' + c, 'w') as f:
        f.write('\n'.join(text))
