__author__ = 'Wenzhi Mao'
__version__ = '1.0.0'

release = [int(x) for x in __version__.split('.')]
del x
__all__ = []


def _StartupPath():
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


def _Startup():
    global _hostname
    import os
    _hostname = os.popen('hostname').read().replace('\n', '')
    if not os.path.exists(os.path.join(_path__, '.Info')):
        os.mkdir(os.path.join(_path__, '.Info'))
    if not os.path.exists(os.path.join(_path__, '.Info', 'com_name.bak')):
        f = open(os.path.join(_path__, '.Info', 'com_name.bak'), 'w')
        f.write(_hostname)
        f.close()
        _clear()
    else:
        f = open(os.path.join(_path__, '.Info', 'com_name.bak'), 'r')
        if f.read() != _hostname:
            f.close()
            f = open(os.path.join(_path__, '.Info', 'com_name.bak'), 'w')
            f.write(_hostname)
            _clear()
        f.close()


def _ABSpath():
    import inspect
    import os.path
    caller_file = inspect.stack()[1][1]
    return os.path.abspath(os.path.dirname(caller_file))


def _make(p):
    from os import path
    from os import popen
    abp = path.abspath(p)
    if path.splitext(path.split(abp)[1])[1] == '.c':
        sop = abp[:-2]+'_c.so'
        popen('mpicc -shared -fPIC -O3 -o '+sop+' '+abp)


_StartupPath()

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
    import os
    dirl = os.listdir(searchpath)
    for d in dirl:
        if os.path.isdir(os.path.join(searchpath, d)):
            _clearSo(os.path.join(searchpath, d))
        elif os.path.isfile(os.path.join(searchpath, d)) and (d.endswith('.so') or d.endswith('pyc')):
            os.remove(os.path.join(searchpath, d))


def _clearData():
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
    _clearSo()
    _clearData()


_Startup()
