# -*- coding: utf-8 -*-
"""This module contains some screen output functions.
"""

from atexit import register as _register


__author__ = 'Wenzhi Mao'

__all__ = ["setVerbo", 'printError', 'printInfo', 'printUpdateInfo',
           'finishUpdate']

_verbo = True

_updating = False


def setVerbo(mark, **kwargs):
    """Set the output on/off."""
    global _verbo
    if mark:
        _verbo = True
    else:
        _verbo = False


def checkIsSublime(**kwargs):
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


def printError(x, **kwargs):
    """Print red error information."""

    global _verbo, _updating
    if _verbo:
        if checkIsSublime():
            print '* {0}'.format(str(x))
        else:
            if _updating:
                finishUpdate()
                _updating = False
            print '\x1b[1;31m* {0}\x1b[0;29m'.format(str(x))


def printInfo(x, **kwargs):
    """Print normal information."""

    global _verbo, _updating
    if _verbo:
        if checkIsSublime():
            print '* {0}'.format(str(x))
        else:
            if _updating:
                finishUpdate()
                _updating = False
            print '\x1b[0;29m* {0}\x1b[0;29m'.format(str(x))


def printUpdateInfo(x, **kwargs):

    global _verbo, _updating
    if _verbo:
        if checkIsSublime():
            print '* {0}'.format(str(x))
        else:
            from os import popen
            import sys
            try:
                l = int(os.popen('stty size').read().split()[-1])
            except:
                l = 0
            sys.stdout.write(
                "\r" + " " * l + '\r\x1b[0;29m* {0}\x1b[0;29m'.format(str(x)))
            sys.stdout.flush()
            _updating = True


def finishUpdate(**kwargs):

    global _verbo, _updating
    if _verbo:
        if not checkIsSublime():
            print ''
            _updating = False


def _output_exit():
    global _updating
    if _verbo:
        if not checkIsSublime():
            if _updating:
                finishUpdate()
                _updating = False

_register(_output_exit)
