'''Some sort functions.
'''

__author__ = 'Wenzhi Mao'

__all__ = ['QuickSort']


def _Startup():
    '''Get _path__.'''
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


def QuickSort(l):
    '''This is a quick sort function based on language C.
    Given a list of integer or float numbers, the function will return a sorted list.'''
    from .c_sort import quicksort
    from numpy import array
    l = quicksort(array(l, dtype=float))
    return l

_Startup()
