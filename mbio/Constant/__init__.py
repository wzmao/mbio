# -*- coding: utf-8 -*-
"""This module contains constants in mbio.

	* :func:`.getconstantfunc` - get constant from file by name

"""

__author__ = 'Wenzhi Mao'

__all__ = []


def getconstantfunc(name, **kwargs):
    """Get constants from file by name."""

    from . import __path__ as path
    from numpy import fromfile
    from os.path import join
    from os import listdir

    path = path[0]
    if not name in listdir(path):
        from ..IO.output import printError
        printError("File {0} not exists.".format(name))
        raise ValueError("File {0} not exists.".format(name))
    temp = fromfile(join(path, name))

    return temp
