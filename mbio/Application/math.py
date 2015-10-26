# -*- coding: utf-8 -*-
"""This module contains some math functions.
In the future plan: Eigenvalue, Inverse, Matrix Multiplication,
                    SVD, PCA
"""

__author__ = 'Wenzhi Mao'

__all__ = ['isSquare', 'ANOVA']


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
    """This is a function to inverse a symetric postive definite matrix."""

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


class ANOVA(object):

    """It is a class for ANOVA analysis. Given the analysis data, 
    output the test result.
    1D data supported now. More dimension could be achieved in the future.
    `data` should be n*2 numpy array or list. The first column is the value 
    and the second column is the label."""

    def __init__(self, data=None, **kwargs):
        """Calculate the ANOVA for the data."""

        from ..IO.output import printInfo
        self.data = data
        self.result = None
        self.pvalue = self.f0 = self.fd = self.sst = self.sstreat = self.mstreat = self.sse = self.mse = self.n = None
        if type(data) == type(None):
            self._calculated = False
        else:
            self.performCalculation(**kwargs)

    def performCalculation(self, alpha=0.05, outprint=True, **kwargs):
        """Perform the ANOVA calculation for the data."""

        from ..IO.output import printInfo, printError
        from numpy import array
        from scipy.stats import f as F

        self._calculated = None
        self.pvalue = self.f0 = self.fd = self.sst = self.sstreat = self.mstreat = self.sse = self.mse = self.n = None
        try:
            self.data = array(self.data, dtype=float)
        except:
            printError("The data could not be transfered to numpy.array")
        if self.data.ndim != 2:
            printError("ANOVA class could only support 1D data now.")
            return None
        if self.data.shape[1] != 2:
            printError("The data should be 2 column data.")
            return None

        labelset = set()
        for i in self.data[:, 1]:
            if not i in labelset:
                labelset.add(i)
        labelset = list(labelset)
        labelset.sort()
        printInfo("{} label(s) found".format(len(labelset)))
        muall = self.data[:, 0].mean()
        sst = ((self.data[:, 0] - muall)**2).sum()
        n = self.data.shape[0]
        ns = array([(self.data[:, 1] == i).sum()
                    for i in labelset], dtype=float)
        mus = array([self.data[:, 0][
                    (self.data[:, 1] == i)].mean() - muall for i in labelset], dtype=float)
        sstreat = (mus**2).dot(ns)
        mstreat = sstreat * 1.0 / (len(ns) - 1)
        mse = (0.0 + sst - sstreat) * 1.0 / (n - len(ns))
        f0 = mstreat / mse
        self.pvalue = 1. - F.cdf(f0, len(ns) - 1, n - len(ns))
        self.f0 = f0
        self.fd = (len(ns) - 1, n - len(ns))
        self.sst = sst
        self.sstreat = sstreat
        self.mstreat = mstreat
        self.sse = (0.0 + sst - sstreat)
        self.mse = mse
        self.n = n
        self._calculated = True
        if outprint:
            printInfo("SS_Total     = {0:13.8f} for {1} data".format(sst, n))
            printInfo("MS_Treatment = {0:13.8f} with {1:6d} of free degrees".format(
                mstreat, self.fd[0]))
            printInfo(
                "MS_Error     = {0:13.8f} with {1:6d} of free degrees".format(mse, self.fd[1]))
            printInfo("F0 = MS_Treatment/MS_Error = {0:12.8f}".format(f0))
            printInfo(
                "p-value      = {0:13.8f} = {1:8.6f}%".format(self.pvalue, self.pvalue * 100))
            if self.pvalue < alpha:
                printInfo(
                    "Reject the null hypothesis at alpha = {}, each class are different.".format(alpha))
            else:
                printInfo(
                    "Accept the null hypothesis at alpha = {}, each class are the same.".format(alpha))
        return None
