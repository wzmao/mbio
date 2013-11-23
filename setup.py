import os
import sys
import platform
from os import sep as dirsep
from os.path import isfile, join

from distutils.core import setup
from distutils.extension import Extension
from distutils.command.install import install

if sys.version_info[:2] < (2, 6):
    sys.stderr.write('Python 2.5 and older is not supported\n')
    sys.exit()

if os.name == 'java':
    sys.stderr.write('JavaOS is not supported\n')
    sys.exit()

try:
    import numpy
except ImportError:
    sys.stderr.write('numpy is not installed, you can find it at: '
                     'http://numpy.scipy.org\n')
    sys.exit()

if [int(dgt) for dgt in numpy.__version__.split('.')[:2]] < [1, 4]:
    sys.stderr.write('numpy v1.4 or later is required, you can find it at: '
                     'http://numpy.scipy.org\n')
    sys.exit()


__version__ = ''
with open('mbio/__init__.py') as inp:
  for line in inp:
      if line.startswith('__version__'):
          exec(line.strip())
          break

PACKAGES = ['mbio',
            'mbio.IO',
            'mbio.Application',
            'mbio.Sequence', ]
PACKAGE_DATA = {
    'mbio': ['Scripts/*.c',
             'Scripts/*.job',]
}

PACKAGE_DIR = {}
for pkg in PACKAGES:
    PACKAGE_DIR[pkg] = join(*pkg.split('.'))

EXTENSIONS = [
    Extension('mbio.Sequence.c_correlation',
              [join('mbio', 'Sequence', 'c_correlation.c'),],
              include_dirs=[numpy.get_include()]),
    Extension('mbio.Sequence.c_shuffle',
              [join('mbio', 'Sequence', 'c_shuffle.c'),],
              include_dirs=[numpy.get_include()]),
    Extension('mbio.Application.c_sort',
              [join('mbio', 'Application', 'c_sort.c'),],
              include_dirs=[numpy.get_include()]),
]

setup(
    name='mbio',
    version=__version__,
    author='Wenzhi Mao',
    author_email='mao.doudoudou@gmail.com',
    description='A Python Package for Biology Self Usage',
    packages=PACKAGES,
    package_dir=PACKAGE_DIR,
    package_data=PACKAGE_DATA,
    ext_modules=EXTENSIONS,
    requires=['NumPy (>=1.5)', ],
    provides=['mbio ({0:s})'.format(__version__)]
)
