__author__ = 'Wenzhi Mao'

__all__ = []


def _Startup():
    from os import path
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()
    Clist = []
    for c in Clist:
        if not path.exists(_path__+'/'+c.replace('.c', '_c.so')):
            from mbio import _make
            _make(_path__+'/'+c)

_Startup()


from . import Fasta
from .Fasta import *
__all__.extend(Fasta.__all__)
