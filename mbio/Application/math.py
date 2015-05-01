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


def eigh(x):
    """This is a function to calculate eigenvalues and eigenvectors."""

    try:
        from scipy.linalg.lapack import dsyevr
        return dsyevr(x)[:2]
    except:
        from numpy.linalg import eigh as n_eigh
        return n_eigh(x)

def invsp(x):
    """This is a function to inverse a sumetric postive definite matrix."""

    try:
        from numpy.linalg import inv
        return inv(x)
    except:
        try:
            from scipy.linalg.lapack import dgetrf, dgetri
            d, e = dgetrf(x)[:2]
            return dgetri(d, e)[0]
        except:
            from ..IO.output import printError
            printError("There is no `inv` function found.")
            return None
