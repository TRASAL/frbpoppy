import glob
import os
import subprocess
import sys

assert sys.version_info >= (3, 0), 'Please run with python3 or higher'


def loc(f):
    """Returns full location path to file"""
    return os.path.join(os.path.dirname(__file__), f)


def run(command):
    """Run a console command given as a list"""
    p = subprocess.Popen(command, stdout=subprocess.PIPE)
    out = p.stdout.read().decode()
    if out:
        print(out)


# Magic to compile fortran code to something called a shared library

all_fortran = glob.glob(loc('./data/models/dm/*.f'))

for f in all_fortran:

    folder = '/'.join(f.split('/')[:-1]) + '/'
    fortran = f.split('/')[-1].split('.')[0]

    setup = ['gfortran',
             '-O2',
             '-fPIC',
             '-fno-second-underscore',
             '-c',
             '-I.',
             '-std=legacy',
             f,
             '-o',
             folder + fortran + '.o',
             ]

    #run(setup)

# Convert .o file to something with which python can interact
gf = ['gfortran',
      '-shared',
      '-o',
      loc('./data/models/dm/libne2001.so'),
      '-fno-second-underscore',
      loc('./data/models/dm/dm.o'),
      loc('./data/models/dm/ne2001.o'),
      loc('./data/models/dm/psr_ne.o'),
      loc('./data/models/dm/dist.o'),
      loc('./data/models/dm/calc_xyz.o'),
      loc('./data/models/dm/density.o'),
      loc('./data/models/dm/glun.o'),
      ]

run(gf)
