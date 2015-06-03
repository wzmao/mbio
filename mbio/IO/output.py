# -*- coding: utf-8 -*-
"""This module contains some screen output functions.
"""

__author__ = 'Wenzhi Mao'

__all__ = ["setVerbo", 'printError', 'printInfo']

_verbo = True


def setVerbo(mark, **kwarg):
    """Set the output on/off."""
    global _verbo
    if mark:
        _verbo = True
    else:
        _verbo = False


def checkIsSublime(**kwarg):
    """Check if the current job is called from sublime."""
    from os import getppid, popen
    try:
        pid = str(getppid())
        f = open('/proc/' + str(pid) + '/cmdline', 'r')
        temp = f.read().replace('\0', ' ')
        f.close()
        if temp.lower().find('sublime') != -1:
            return True
        else:
            return False
    except:
        return False


def printError(x, **kwarg):
    """Print red error information."""

    if _verbo:
        if checkIsSublime():
            print '* {0}'.format(str(x))
        else:
            print '\x1b[1;31m* {0}\x1b[0;29m'.format(str(x))


def printInfo(x, **kwarg):
    """Print normal information."""

    if _verbo:
        if checkIsSublime():
            print '* {0}'.format(str(x))
        else:
            print '\x1b[0;29m* {0}\x1b[0;29m'.format(str(x))
