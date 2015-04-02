# -*- coding: utf-8 -*-
"""This module contains features for IO in mbio.

Screen output
========

	* :func:`.printError` - print red error information
	* :func:`.printInfo` - print normal information

File IO
In the future plan: mat(MATLAB), dcd
========

"""

__author__ = 'Wenzhi Mao'

__all__ = []


from . import output
from .output import *
__all__.extend(output.__all__)
