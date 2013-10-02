__author__ = 'Wenzhi Mao'

__all__ = []


def _Startup():
    '''Get _path__ and compile files.'''
    from os import path
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()
    Clist = ['sort.c']
    Compiler = ['']
    Options = ['']
    Libs = ['']
    for c, cp, o, l in zip(Clist, Compiler, Options, Libs):
        if not path.exists(path.join(_path__, c.replace('.c', '_c.so'))):
            from mbio.Application import compile
            compile.make(path.join(_path__, c), cp, o, l)

_Startup()


from . import sort
from .sort import *
__all__.extend(sort.__all__)

from . import cluster
# from .cluster import *
# __all__.extend(cluster.__all__)

from . import job_organization
from .job_organization import *
__all__.extend(job_organization.__all__)

from . import math
# from .math import *
# __all__.extend(math.__all__)