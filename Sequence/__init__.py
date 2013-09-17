__author__ = 'Wenzhi Mao'

__all__ = []


def _Startup():
    from os import path
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()
    Clist = ['correlation.c']
    for c in Clist:
        if not path.exists(path.join(_path__, c.replace('.c', '_c.so'))):
            from mbio import _make
            _make(path.join(_path__, c))

_Startup()


from . import correlation
from .correlation import *
__all__.extend(correlation.__all__)

from . import shuffle
from .shuffle import *
__all__.extend(shuffle.__all__)
