__author__ = 'Wenzhi Mao'

__all__ = []


def _Startup():
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()
    from os import path
    Clist = ['mi.c', 'omes.c']
    for c in Clist:
        if not path.exists(_path__+'/'+c.replace('.c', '_c.so')):
            from mbio import _make
            _make(_path__+'/'+c)

_Startup()

from . import MI
from .MI import *
__all__.extend(MI.__all__)

from . import OMES
from .OMES import *
__all__.extend(OMES.__all__)
