__author__ = 'Wenzhi Mao'
__version__ = '1.0.0'

release = [int(x) for x in __version__.split('.')]
del x
__all__ = ['_make', '_ABSpath']


def _Startup():
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()
    import os
    _hostname = os.popen('hostname').read().replace('\n','')
    if not os.path.exists(_path__+'/.Info'):
        os.mkdir(_path__+'/.Info')
    if not os.path.exists(_path__+'/.Info/com_name'):
        f = open(_path__+'/.Info/com_name','w')
        f.write(_hostname)
        f.close()
        _clearSo(_path__)
        _clearData()
    else:
        f = open(_path__+'/.Info/com_name','r')
        if f.read() != _hostname:
            f.close()
            f = open(_path__+'/.Info/com_name','w')
            f.write(_hostname)
            _clearSo(_path__)
            _clearData()
        f.close()


def _make(p):
    from os import path
    from os import popen
    abp = path.abspath(p)
    if path.splitext(path.split(abp)[1])[1] == '.c':
        sop = abp[:-2]+'_c.so'
        popen('mpicc -shared -fPIC -O3 -o '+sop+' '+abp)


def _ABSpath():
    import inspect
    import os.path
    caller_file = inspect.stack()[1][1]
    return os.path.abspath(os.path.dirname(caller_file))


def _clearSo(searchpath):
    import os
    dirl = os.listdir(searchpath)
    for d in dirl:
        if os.path.isdir(os.path.join(searchpath, d)):
            _clearSo(os.path.join(searchpath, d))
        elif os.path.isfile(os.path.join(searchpath, d)) and d.endswith('.so'):
            os.remove(os.path.join(searchpath, d))


def _clearData():
    import os
    if os.path.exists(_path__+'/.Cache/'):
        for i in os.listdir(_path__+'/.Cache/'):
            os.remove(_path__+'/.Cache/'+i)
        os.removedirs(_path__+'/.Cache/')
    if os.path.exists(_path__+'/.Result/'):
        for i in os.listdir(_path__+'/.Result/'):
            os.remove(_path__+'/.Result/'+i)
        os.removedirs(_path__+'/.Result/')


_Startup()

from . import Sort
from .Sort import *
__all__.extend(Sort.__all__)
__all__.append('Sort')


from . import Correlation
from .Correlation import *
__all__.extend(Correlation.__all__)
__all__.append('Correlation')


from . import IO
from .IO import *
__all__.extend(IO.__all__)
__all__.append('IO')