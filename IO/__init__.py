__author__ = 'Wenzhi Mao'

__all__ = []


def _Startup():
    '''Get _path__ and compile files.'''
    from os import path
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()
    Clist = []
    Compiler = []
    Options = []
    for c, cp, o in zip(Clist,Compiler,Options):
        if not path.exists(path.join(_path__, c.replace('.c', '_c.so'))):
            from mbio.Application import compile
            compile.make(path.join(_path__, c), cp, o)

_Startup()


from . import fasta
from .fasta import *
__all__.extend(fasta.__all__)

from . import matrix
from .matrix import *
__all__.extend(matrix.__all__)

from . import error
from .error import *
__all__.extend(error.__all__)