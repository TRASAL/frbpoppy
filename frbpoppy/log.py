"""Quick and dirty logging functions."""
import inspect


def pprint(*s, output=True):
    """Hack to make for more informative print statements."""
    f = inspect.stack()[1][1].split('/')[-1]
    m = '{:13.13} |'.format(f)

    if output:
        print(m, *s)
    else:
        lines = []
        for e in s:
            lines.append('\n'.join([f'{m} {f}' for f in e.split('\n')]))
        return '\n'.join(lines)
