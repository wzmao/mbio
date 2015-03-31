__author__ = 'Wenzhi Mao'
__version__ = '1.0.0'

release = [int(x) for x in __version__.split('.')]
del x, release
__all__ = []


def _ABSpath():
    '''Get absolute path for the script.'''
    import inspect
    import os.path
    caller_file = inspect.stack()[1][1]
    return os.path.abspath(os.path.dirname(caller_file))


from . import Application
from .Application import *
__all__.extend(Application.__all__)
__all__.append('Application')


from . import Sequence
from .Sequence import *
__all__.extend(Sequence.__all__)
__all__.append('Sequence')


from . import IO
from .IO import *
__all__.extend(IO.__all__)
__all__.append('IO')


def _clearSo(searchpath=_path__):
    '''Delete all .so file and .pyc in mbio.
    Could affect some functions. Restart the program for normal usage.'''
    import os
    dirl = os.listdir(searchpath)
    for d in dirl:
        if os.path.isdir(os.path.join(searchpath, d)):
            _clearSo(os.path.join(searchpath, d))
        elif os.path.isfile(os.path.join(searchpath, d)) and (d.endswith('.so') or d.endswith('pyc')):
            os.remove(os.path.join(searchpath, d))


def _clearData():
    '''Delete all running job files.'''
    import os
    if os.path.exists(os.path.join(_path__, '.Cache')):
        for i in os.listdir(os.path.join(_path__, '.Cache')):
            os.remove(os.path.join(_path__, '.Cache', i))
        os.removedirs(os.path.join(_path__, '.Cache'))
    if os.path.exists(os.path.join(_path__, '.Result')):
        for i in os.listdir(os.path.join(_path__, '.Result')):
            os.remove(os.path.join(_path__, '.Result', i))
        os.removedirs(os.path.join(_path__, '.Result'))


def _clear():
    '''Make the package back to start status.'''
    _clearSo()
    _clearData()
