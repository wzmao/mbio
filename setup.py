import os
import sys
import platform
from os import sep as dirsep
from os.path import isfile, join

from distutils.core import setup
from distutils.extension import Extension
from distutils.command.install import install

# os.environ['CC']='gcc'
# os.environ['LDSHARED']='g++'
# os.environ['CXX']='g++'
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

if [int(dgt) for dgt in numpy.__version__.split('.')[:2]] < [1, 9]:
    sys.stderr.write('numpy v1.9 or later is required, you can find it at: '
                     'http://numpy.scipy.org\n')
    sys.exit()

try:
    import scipy
except ImportError:
    sys.stderr.write('scipy is not installed, you can find it at: '
                     'http://www.scipy.org/\n')
    sys.exit()

if [int(dgt) for dgt in scipy.__version__.split('.')[:2]] < [0, 15]:
    sys.stderr.write('scipy v0.15 or later is required, you can find it at: '
                     'http://numpy.scipy.org\n')
    sys.exit()

__version__ = ''
with open('mbio/__init__.py') as inp:
    for line in inp:
        if line.startswith('__version__'):
            exec(line.strip())
            break

PACKAGES = ['mbio',
            'mbio.Application',
            'mbio.Constant',
            'mbio.EM',
            'mbio.IO',
            'mbio.Sequence', ]

PACKAGE_DATA = {'mbio.Constant': ['pplus', 'pminus', 'psigplus', 'psigminus']}

PACKAGE_DIR = {}
for pkg in PACKAGES:
    PACKAGE_DIR[pkg] = join(*pkg.split('.'))

optimize = ['-O3'] if platform.system() == 'Linux' else []
zlib = ['zdll.lib'] if platform.system() == 'Windows' else []

EXTENSIONS = [
    Extension('mbio.Sequence.Ccorrelation',
              [join('mbio', 'Sequence', 'Ccorrelation.c'), ],
              include_dirs=[numpy.get_include()],
              extra_compile_args=[] + optimize,
              extra_link_args=[] + optimize),
    Extension('mbio.Sequence.Ccorrelation_p',
              [join('mbio', 'Sequence', 'Ccorrelation_p.c'), ],
              include_dirs=[numpy.get_include()],
              extra_compile_args=['-fopenmp'] + optimize,
              extra_link_args=['-lgomp'] + optimize),
    Extension('mbio.Sequence.Cshuffle',
              [join('mbio', 'Sequence', 'Cshuffle.c'), ],
              include_dirs=[numpy.get_include()],
              extra_compile_args=[] + optimize,
              extra_link_args=[] + optimize),
    Extension('mbio.Application.c_algorithm',
              [join('mbio', 'Application', 'c_algorithm.c'), ],
              include_dirs=[numpy.get_include()],
              extra_compile_args=[] + optimize,
              extra_link_args=[] + optimize),
    Extension('mbio.Sequence.Cfasta',
              [join('mbio', 'Sequence', 'Cfasta.c'), ],
              include_dirs=[numpy.get_include()],
              extra_compile_args=[] + optimize,
              extra_link_args=[] + optimize),
    Extension('mbio.EM.Cmrc',
              [join('mbio', 'EM', 'Cmrc.c'), ],
              include_dirs=[numpy.get_include()],
              extra_compile_args=[] + optimize,
              extra_link_args=[] + optimize + zlib),
    Extension('mbio.EM.Cmrcmodel_p',
              [join('mbio', 'EM', 'Cmrcmodel_p.c'), ],
              include_dirs=[numpy.get_include()],
              extra_compile_args=['-fopenmp'] + optimize,
              extra_link_args=['-lfftw3f', '-lgomp'] + optimize),
    # Extension('mbio.Learning.CNN_p',
    #           [join('mbio', 'Learning', 'CNN_p.c'), ],
    #           include_dirs=[numpy.get_include()],
    #           extra_compile_args=['-fopenmp', '-O3'],
    #           extra_link_args=['-lgomp', '-O3']),
]

setup(
    name='mbio',
    version=__version__,
    author='Wenzhi Mao',
    author_email='maowenzhi@yandex.com',
    description='A Python Package for Biology from Wenzhi Mao',
    long_description='''A Python Package for Biology from Wenzhi Mao
    Include Sequence analysis tools and MRC tools.''',
    packages=PACKAGES,
    package_dir=PACKAGE_DIR,
    package_data=PACKAGE_DATA,
    ext_modules=EXTENSIONS,
    requires=['NumPy (>=1.9)', 'Scipy (>=0.15)'],
    provides=['mbio ({0:s})'.format(__version__)],
    url='https://github.com/wzmao/mbio',
    license='MIT',
    platforms=['Linux', 'Windows'],
)
