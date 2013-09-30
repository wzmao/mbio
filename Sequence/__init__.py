__author__ = 'Wenzhi Mao'

__all__ = []


def _Startup():
    '''Get _path__ and compile files.'''
    from os import path
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()
    Clist = ['correlation.c']
    Compiler = ['']
    Options = ['']
    for c, cp, o in zip(Clist,Compiler,Options):
        if not path.exists(path.join(_path__, c.replace('.c', '_c.so'))):
            from mbio.Application import compile
            compile.make(path.join(_path__, c), cp, o)

_Startup()


from . import correlation
from .correlation import *
__all__.extend(correlation.__all__)

from . import shuffle
from .shuffle import *
__all__.extend(shuffle.__all__)
