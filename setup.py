"""Script to do a local install of frbpoppy."""

import os
import sys
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

        all_fortran = glob.glob(loc('./data/models/ne2001/*.f'))

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
              loc('./data/models/ne2001/libne2001.so'),
              '-fno-second-underscore',
              loc('./data/models/ne2001/dm.o'),
              loc('./data/models/ne2001/ne2001.o'),
              loc('./data/models/ne2001/psr_ne.o'),
              loc('./data/models/ne2001/dist.o'),
              loc('./data/models/ne2001/calc_xyz.o'),
              loc('./data/models/ne2001/density.o'),
              loc('./data/models/ne2001/glun.o'),
              ]

        check_call(gf)

        # Some additional information of users of Python 3.5 or older
        v = sys.version_info
        if v[0] == 3 and v[1] < 6:
            m = 'Note from frbpoppy\n'
            m += '==================\n'
            m += ' - It seems you are using an old Python version (v{}.{}) \n'
            m = m.format(v[0], v[1])
            m += (
                 ' - We recommend upgrading to >=3.6 if possible \n'
                 " - As a work around we have prepared a bash script to "
                 'backport frbpoppy to versions 3.0<=v<=3.5 \n'
                 ' - This will pip3 install a library called future_fstrings '
                 'and remove all f-strings from the code '
                 )
            print(m)
            i = input(' > Would you like frbpoppy to run this script? [y/n]')
            if i == 'y' or i == 'Y':
                c = ['bash', 'backport.sh']
                check_call(c)
            else:
                ' - Avoiding runing bash script.'

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
      python_requires='>=3.4',
      install_requires=['bokeh >= 1.3.4',
                        'numpy',
                        'pandas >= 0.23.4',
                        'scipy >= 1.1.0',
                        'SQLAlchemy >= 1.3.0',
                        'matplotlib >= 2.2.3,<3.1',
                        'requests >= 2.20.0.',
                        'dill >= 0.3.1.1',
                        'tqdm',
                        'joblib'],
      cmdclass={'develop': PostDevelopCommand})
