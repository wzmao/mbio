# -*- coding: utf-8 -*-
"""This module contains some screen output functions.
"""

__author__ = 'Wenzhi Mao'

__all__ = ['printError', 'printInfo']


def printError(x):
    """Print red error information."""

    print '\x1b[1;31m* {0}\x1b[0;29m'.format(str(x))


def printInfo(x):
    """Print normal information."""

    print '\x1b[0;29m* {0}\x1b[0;29m'.format(str(x))
