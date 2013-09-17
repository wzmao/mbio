__author__ = 'Wenzhi Mao'

__all__ = []


def _Startup():
    from os import path
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()
    Clist = ['sort.c']
    for c in Clist:
        if not path.exists(path.join(_path__, c.replace('.c', '_c.so'))):
            from mbio import _make
            _make(path.join(_path__, c))

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