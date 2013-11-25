__author__ = 'Wenzhi Mao'

__all__ = []


def _Startup():
    '''Get _path__ and compile files.'''
    from os import path
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


_Startup()


from . import correlation
from .correlation import *
__all__.extend(correlation.__all__)

from . import shuffle
from .shuffle import *
__all__.extend(shuffle.__all__)

from . import msa
from .msa import *
__all__.extend(msa.__all__)

from . import msaio
from .msaio import *
__all__.extend(msaio.__all__)

from . import analysis
from .analysis import *
__all__.extend(analysis.__all__)
