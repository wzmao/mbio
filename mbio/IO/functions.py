# -*- coding: utf-8 -*-
"""This module contains some general IO functions.
"""

__author__ = 'Wenzhi Mao'

__all__ = ["parseMRC", "parseMRCHeader", "writeMRC"]


def parseMRC(filename=None, **kwargs):
    """Parse the MRC from a file."""

    from .mrc import MRC
    from os.path import exists, isfile

    if type(filename) == type(None):
        from .output import printError
        printError("The filename is wrong.")
        return None
    if exists(filename) and isfile(filename):
        return MRC(filename=filename)
    else:
        from .output import printError
        printError("The filename doesn't exists or is not a file.")


def parseMRCHeader(filename=None, **kwargs):
    """Parse the MRC header from a file."""
    from .mrc import MRCHeader

    from os.path import exists, isfile

    if type(filename) == type(None):
        from .output import printError
        printError("The filename is wrong.")
        return None
    if exists(filename) and isfile(filename):
        return MRCHeader(filename=filename)
    else:
        from .output import printError
        printError("The filename doesn't exists or is not a file.")


def writeMRC(filename, mrc, **kwargs):
    """Write the MRC file into a file."""
    from .mrc import MRC
    if not isinstance(mrc, MRC):
        from .output import printError
        printError("The mrc is not a mbio MRC instance.")
        return None
    mrc.writeData(filename=filename, **kwargs)
