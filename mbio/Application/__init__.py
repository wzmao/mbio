__author__ = 'Wenzhi Mao'

__all__ = []


def _Startup():
    '''Get _path__ and compile files.'''
    from os import path
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


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
