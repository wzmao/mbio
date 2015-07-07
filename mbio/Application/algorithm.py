# -*- coding: utf-8 -*-
"""This module contains some algorithm functions.
"""

__author__ = 'Wenzhi Mao'

__all__ = ['quickSort', 'binarySearch']


def quickSort(l):
    """This is a quick sort function based on language C.
    Given a list of integer or float numbers, the function will return a sorted array.
    It will not change the given variable directly.
    It has faster speed than the default python sort function."""

    from .c_algorithm import quicksort
    from numpy import array
    l = quicksort(array(l, dtype=float))
    return l


def binarySearch(se, v, s=None, e=None):
    """This is a binary search algorithm.
    Now it is in python, we will change it to C some time later.
    In a sorted list, find the first number index >= the target."""

    if type(s) == type(None):
        s = 0
    if type(e) == type(None):
        e = len(se) - 1
    if se[s] >= v:
        return s
    if e - s == 0:
        return s
    if se[(s + e) // 2] == v:
        return (s + e) // 2
    if e - s == 1:
        return e
    if se[(s + e) // 2] < v:
        return binarySearch(se, v, (s + e) // 2 + 1, e)
    if se[(s + e) // 2] > v:
        return binarySearch(se, v, s, (s + e) // 2)
