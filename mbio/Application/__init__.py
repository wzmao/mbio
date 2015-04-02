"""This module contains features for general applications in mbio.

Math (If it grows too big in the future, make it alone)
========

	* :func:`.isSquare` - check if a integer is a square

Sort
========

	* :func:`.quickSort` - Quick Sort a list of float/int

"""

__author__ = 'Wenzhi Mao'

__all__ = []


from . import sort
from .sort import *
__all__.extend(sort.__all__)

# from . import cluster
# from .cluster import *
# __all__.extend(cluster.__all__)

# from . import job_organization
# from .job_organization import *
# __all__.extend(job_organization.__all__)

from . import math
# from .math import *
# __all__.extend(math.__all__)
