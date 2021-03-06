# -*- coding: utf-8 -*-
"""This module contains features for general applications in mbio.

Math (If it grows too big in the future, make it alone)
========

    * :func:`.isSquare` - check if a integer is a square
    * :func:`.eigh` - return the eigenvalues and eigenvectors
    * :func:`.invsp` - inverse a symetric postive definite matrix
    * :func:`.performRegression - perform the regression and t-test for given data
    * :func:`.performPolyRegression - perform polynomial regression and F-test for given data

    * :class:`.ANOVA` - a class to perform the ANOVA analysis

Algorithm
========

    * :func:`.quickSort` - Quick Sort a list of float/int
    * :func:`.binarySearch` - binary search in a sorted list-like object

"""

__author__ = 'Wenzhi Mao'

__all__ = []


from . import algorithm
from .algorithm import *
__all__.extend(algorithm.__all__)

# from . import cluster
# from .cluster import *
# __all__.extend(cluster.__all__)

# from . import job_organization
# from .job_organization import *
# __all__.extend(job_organization.__all__)

from . import math
from .math import *
__all__.extend(math.__all__)

from . import setting
from .setting import *
__all__.extend(setting.__all__)

from . import plotting
from .plotting import *
# __all__.extend(plotting.__all__)
