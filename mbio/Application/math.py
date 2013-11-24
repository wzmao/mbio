'''Some math functions.
'''

__author__ = 'Wenzhi Mao'

__all__ = ['issquare']


def _Startup():
    '''Get _path__.'''
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


def issquare(x):
    '''It is a function to determine if the integer is a square of another integer.'''
    try:
        xi = int(x)
    except:
        return None
    if xi != x:
        from mbio.IO import error
        error.errorprint('The number is not integer.')
        return None
    if x < 0:
        from mbio.IO import error
        error.errorprint('The number is negative.')
        return None
    x = xi
    sq = x ** .5
    if abs(int(round(sq, 0)) ** 2 - x) < 1e-10:
        return True
    else:
        return False


_Startup()
