'''It is a module to compile C files.'''


__author__ = 'Wenzhi Mao'
__all__ = []


def _Startup():
    '''Get _path__.'''
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


def make(p,):
    '''Compile C files using mpicc.'''
    from os import path
    from os import popen
    abp = path.abspath(p)
    if path.splitext(path.split(abp)[1])[1] == '.c':
        sop = abp[:-2]+'_c.so'
        popen('gcc -shared -fPIC -O3 -o '+sop+' '+abp)