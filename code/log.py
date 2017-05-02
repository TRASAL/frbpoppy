import inspect


def pprint(*s):
    """Hack to make for more informative print statements"""
    f = inspect.stack()[1][1].split('/')[-1]
    m = '{:13.13} |'.format(f)
    print(m, *s)

