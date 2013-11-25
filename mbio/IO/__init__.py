__author__ = 'Wenzhi Mao'

__all__ = []


def _Startup():
    '''Get _path__ and compile files.'''
    from os import path
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


_Startup()


from . import matrix
from .matrix import *
__all__.extend(matrix.__all__)

from . import error
from .error import *
__all__.extend(error.__all__)
