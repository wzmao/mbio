# -*- coding: utf-8 -*-
"""This module contains features for EM and MRC in mbio.

Classes
=======

    * :class:`.MRC` - store MRC header and data and some related functions
    * :class:`.MRCHeader` - store MRC header

EM Analysis functions
=======

    * :func:`.interpolationball` - interpolate at a position as a ball in EM map
    * :func:`.interpolationcube` - interpolate at a position as a cube in EM map
    * :func:`.genPvalue` - generate p-value for a structure/coordinates in EM data
    * :func:`.calcPcutoff` - calculate the p-value cutoff bu recognize specific pattern
    * :func:`.showPcutoff` - plot the p-value cutoff bu recognize specific pattern
    * :func:`.genPvalueSample` - generate p-value sample for structures and EM

"""

__author__ = 'Wenzhi Mao'

__all__ = []

from . import analysis
from .analysis import *
__all__.extend(analysis.__all__)

from . import mrc
from .mrc import *
__all__.extend(mrc.__all__)
