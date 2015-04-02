# -*- coding: utf-8 -*-
"""This module contains some math functions.
In the future plan: Eigenvalue, Inverse, Matrix Multiplication,
                    SVD, PCA
"""

__author__ = 'Wenzhi Mao'

__all__ = ['isSquare']


def isSquare(x):
    """It is a function to determine if the given integer is a square integer."""

    try:
        xi = int(x)
    except:
        return None
    if xi != x:
        from ..IO.output import printError
        printError('The number is not integer.')
        return None
    if x < 0:
        from ..IO.output import printError
        printError('The number is negative.')
        return None
    x = xi
    sq = x ** .5
    if abs(int(round(sq, 0)) ** 2 - x) < 1e-10:
        return True
    else:
        return False
