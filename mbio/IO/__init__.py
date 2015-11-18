# -*- coding: utf-8 -*-
"""This module contains features for IO in mbio.

Screen output
========

	* :func:`.printError` - print red error information
	* :func:`.printInfo` - print normal information
	* :func:`.setVerbo` - set the printout on/off

File IO
In the future plan: mat(MATLAB), dcd
========

  MRC Realted
	* :func:`.parseMRC` - parse MRC header and data from a MRC file
	* :func:`.parseMRCHeader` - parse MRC header from a MRC file
	* :func:`.writeMRC` - write the MRC into a file
  PDB Related
	* :func:`.parsePDB` - parse PDB header and data from a PDB file
	* :func:`.writePDB` - write the PDB into a file

"""

__author__ = 'Wenzhi Mao'

__all__ = []

from . import output
from .output import *
__all__.extend(output.__all__)

from . import PDBio
from .PDBio import *
__all__.extend(PDBio.__all__)

from . import MRCio
from .MRCio import *
__all__.extend(MRCio.__all__)
