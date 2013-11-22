'''Some error and information print.
'''

__author__ = 'Wenzhi Mao'

__all__ = ['errorprint']


def _Startup():
    '''Get _path__.'''
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


def errorprint(x):
    '''It is a print red error information function.'''
    print '\x1b[1;31m* {0}\x1b[0;29m'.format(str(x))


_Startup()
