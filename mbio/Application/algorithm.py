# -*- coding: utf-8 -*-
"""This module contains some algorithm functions.
"""

__author__ = 'Wenzhi Mao'

__all__ = ['quickSort']


def quickSort(l):
    """This is a quick sort function based on language C.
    Given a list of integer or float numbers, the function will return a sorted array.
    It will not change the given variable directly.
    It has faster speed than the default python sort function."""

    from .c_algorithm import quicksort
    from numpy import array
    l = quicksort(array(l, dtype=float))
    return l
