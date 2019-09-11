"""Script to do a local install of frbpoppy."""

import os
import glob
from setuptools import setup
from setuptools.command.develop import develop
from subprocess import check_call


class PostDevelopCommand(develop):
    """Post-installation for development mode."""

    def run(self):
        """Compile fortran libraries for NE2001."""
        def loc(f):
            """Full location path to file."""
            return os.path.join(os.path.dirname(__file__), f)

        all_fortran = glob.glob(loc('./data/models/dm/*.f'))

        for f in all_fortran:

            folder = '/'.join(f.split('/')[:-1]) + '/'
            fortran = f.split('/')[-1].split('.')[0]

            to_o = ['gfortran',
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

            check_call(to_o)

        if os.name == 'mac':  # Mac
            flag = '-dynamiclib'
        if os.name == 'nt':  # Windows
            flag = '-shared'
        else:  # Linux
            flag = '-shared'

        # Convert .o file to something with which python can interact
        gf = ['gfortran',
              flag,
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

        check_call(gf)

        develop.run(self)


setup(name='frbpoppy',
      version='1.0.0',
      description='Fast Radio Burst Population Synthesis',
      long_description=open('README.rst').read(),
      url='http://github.com/davidgardenier/frbpoppy',
      author='David Gardenier',
      author_email='gardenier@astron.nl',
      license='MIT',
      packages=['frbpoppy'],
      zip_safe=False,
      python_requires='>=3.0',
      install_requires=['bokeh >= 1.3.4',
                        'numpy >= 1.17.0',
                        'pandas >= 0.23.4',
                        'scipy >= 1.1.0',
                        'SQLAlchemy >= 1.3.0',
                        'matplotlib >= 2.2.3',
                        'requests >= 2.20.0.',
                        'future-fstrings >= 1.2.0'],
      cmdclass={'develop': PostDevelopCommand})
