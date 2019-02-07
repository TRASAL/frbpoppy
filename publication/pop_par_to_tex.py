""" Get population parameters, and convert into Latex table."""

from frbpoppy import paths

path = paths.code + '/../tests/quick.py'

with open(path, 'r') as f:
    lines = f.readlines()

# Parse the arguments given to a standard population
arguments = False
for line in lines:
    if 'def standard' in line:
        arguments = True

    if arguments:
        line = line.strip()
        line = line.strip(',')
        line = line.strip(')')

        try:
            par, val = line.split('=')
            print(par, '|', val)
        except ValueError:
            pass

        if 'return' in line:
            arguments = False
