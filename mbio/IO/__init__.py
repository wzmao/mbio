# -*- coding: utf-8 -*-
"""This module contains features for IO in mbio.

Classes
=======

    * :class:`.MRC` - store MRC header and data and some related functions
    * :class:`.MRCHeader` - store MRC header

Screen output
========

	* :func:`.printError` - print red error information
	* :func:`.printInfo` - print normal information

File IO
In the future plan: mat(MATLAB), dcd
========
	* :func:`.parseMRC` - parse MRC header and data from a MRC file
	* :func:`.parseMRCHeader` - parse MRC header from a MRC file
	* :func:`.writeMRC` - write the MRC into a file

"""

__author__ = 'Wenzhi Mao'

__all__ = []

from . import functions
from .functions import *
__all__.extend(functions.__all__)

from . import output
from .output import *
__all__.extend(output.__all__)

from . import mrc
from .mrc import *
__all__.extend(mrc.__all__)
