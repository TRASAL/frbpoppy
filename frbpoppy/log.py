"""Quick and dirty logging functions."""
import inspect
import sys


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


def progressbar(it, prefix="", size=69, file=sys.stdout):
    """Progressbar from adapted from Stack Overflow.

    Args:
        it (generator): range of values
        prefix (str): Words displayed before the progress bar
        size (int): Display width
        file: Where to direct output

    Returns:
        type: Description of returned object.

    """
    count = len(it)
    size -= len(prefix)

    def show(j):
        x = int((size)*j/count)
        print(f'{prefix} [{"#"*x}{"."*(size-x)}] {j}/{count}')

    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)

    file.flush()
