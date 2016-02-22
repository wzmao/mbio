# -*- coding: utf-8 -*-
"""This module contains some math and statistics functions.
In the future plan: Eigenvalue, Inverse, Matrix Multiplication,
                    SVD, PCA
"""

__author__ = 'Wenzhi Mao'

__all__ = ['isSquare', 'ANOVA', 'performRegression', 'performPolyRegression']


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
        sst = ((self.data[:, 0] - muall) ** 2).sum()
        n = self.data.shape[0]
        ns = array([(self.data[:, 1] == i).sum()
                    for i in labelset], dtype=float)
        mus = array([self.data[:, 0][
                    (self.data[:, 1] == i)].mean() - muall for i in labelset], dtype=float)
        sstreat = (mus ** 2).dot(ns)
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


def performRegression(x, y, const=True, alpha=0.05, label=None, **kwargs):
    """Make regression analysis of array data. And test each parameter using t-test.

    `x` must be a N*a array. `y` must be a N*1 array.
    If `x` or `y` just has one dimension, it could be a 1D array and converted automatically.

    `const` is `True` default and it will detect the are there constant in `x`.
    If no constant in `x`, it will add a new column at the end.

    `alpha` is used to test each parameter.

    `label` could be used for output."""

    from numpy import ndarray, array, hstack, ones
    from numpy.linalg.linalg import inv
    from ..IO.output import printError, printInfo
    from scipy.stats import t

    if not isinstance(x, ndarray) or not isinstance(y, ndarray):
        try:
            x = array(x, dtype=float)
            y = array(y, dtype=float)
        except:
            printError(
                "x and y must be numpy array or could be converted to numpy array.")
            return None
    x = array(x, dtype=float)
    y = array(y, dtype=float)
    if x.ndim == 2:
        pass
    elif x.ndim == 1:
        x.resize((x.size, 1))
    else:
        printError("x must be 1D or 2D data.")
        return None
    if y.ndim == 2:
        if y.shape[1] != 1:
            printInfo("Just take the first column of y.")
            y = y[:, 0:1]
    elif y.ndim == 1:
        y.resize((y.size, 1))
    else:
        printError("y must be 1D or 2D data.")
        return None
    if x.shape[0] != y.shape[0]:
        printError("x and y must have same first dimension.")
        return None
    if type(label) == type(None):
        label = ['x' + str(i + 1) for i in xrange(x.shape[1])]
    else:
        label = [str(i) for i in label]
    if len(label) != x.shape[1]:
        printError(
            "The length of label does not match data. Dismiss the label.")
        label = ['x' + str(i + 1) for i in xrange(x.shape[1])]

    addconst = 0
    if const:
        hasconst = False
        for i in xrange(x.shape[1]):
            if len(set(x[:, i])) == 1:
                hasconst = True
                break
        if not hasconst:
            x = hstack((x, ones((x.shape[0], 1))))
            addconst = 1
            label.append('c')
            printInfo(
                "Add const automatically. If you don't want to add const, use `const = False`")

    cov = inv(x.T.dot(x))
    beta = cov.dot(x.T).dot(y)
    r = y - x.dot(beta)
    sigma2 = ((r.T.dot(r)) / (x.shape[0] - x.shape[1]))[0, 0]
    if sigma2 == 0:
        sigma2 = 5e-324
    st = '\ty = '
    for i in xrange(x.shape[1] - 1):
        st += "{0:+10.6f}*{1:s} ".format(beta[i, 0], label[i])
    if addconst:
        st += "{0:+10.6f}".format(beta[-1, 0])
    else:
        st += "{0:+10.6f}*{1:s}".format(beta[-1, 0], label[i + 1])
    printInfo("The result is :")
    printInfo(st)
    printInfo("Test each parameter.")
    printInfo("\t{0:^5s}{1:^15s}{2:^15s}{3:^15s}{4:^5s}{5:^9s}{6:^5s}".format(
        "xi", "Para", "Sigma", "t-statistics", 'FD', "p-value", 'Sig'))
    p = []
    ts = []
    sig = []
    sigma = []
    for i in xrange(x.shape[1]):
        sigma.append((sigma2 * cov[i, i]) ** .5)
        ts.append(beta[i][0] / sigma[-1])
        p.append((1. - t.cdf(abs(ts[-1]), x.shape[0] - x.shape[1])) * 2)
        sig.append("Yes" if 2. * (1. - t.cdf(abs(beta[i][0] / (
            (sigma2 * cov[i, i]) ** .5)), x.shape[0] - x.shape[1])) < alpha else 'No')
        printInfo("\t{0:^5s}{1:^15.6e}{2:^15.6e}{3:^15.6e}{4:^5d}{5:^9f}"
                  "{6:^5s}".format(label[i],
                                   beta[i][0],
                                   sigma[-1],
                                   ts[-1],
                                   x.shape[0] - x.shape[1],
                                   p[-1],
                                   sig[-1]))
    p = array(p)
    ts = array(ts)
    sig = array(sig)
    sigma = array(sigma)
    return {'beta': beta, 'p': p, 't': ts, "label": label, 'sig': sig, 'sigma': sigma}


def performPolyRegression(y, degree=2, **kwargs):
    '''Build regression with higher degree of polynomial.
    Use the index to build the polynomial.

    The *orthogonal unit* scale is used up to 4 degrees. const is not included.
    You could specific the `const=False` to disable the const.
    '''

    from numpy import ndarray, array, arange, zeros
    from ..IO.output import printError, printInfo

    if not isinstance(y, ndarray):
        try:
            y = array(y, dtype=float)
        except:
            printError(
                "y must be numpy array or could be converted to numpy array.")
            return None
    y = array(y, dtype=float)
    if y.ndim == 2:
        if y.shape[1] != 1:
            printInfo("Just take the first column of y.")
            y = y[:, 0:1]
    elif y.ndim == 1:
        y.resize((y.size, 1))
    else:
        printError("y must be 1D or 2D data.")
        return None
    if not degree in [1, 2, 3, 4]:
        printError("degree must between 1 and 4.")
    if degree + 1 >= y.shape[0]:
        printError("The degree must be less than the data size.")
        return None

    k = y.shape[0] * 1.0
    poly = zeros((k, degree))
    t = arange(k, dtype=float)
    t = t - t.mean()
    label = []
    kwargs.pop('label', None)
    for i in xrange(degree):
        if i == 0:
            label.append('x')
        else:
            label.append('x^' + str(i + 1))
        if i == 0:
            poly[:, i] = t
            poly[:, i] = poly[:, i] / ((poly[:, i] ** 2).sum()) ** .5
        elif i == 1:
            poly[:, i] = t ** 2 - (k ** 2. - 1) / 12
            poly[:, i] = poly[:, i] / ((poly[:, i] ** 2).sum()) ** .5
        elif i == 2:
            poly[:, i] = t ** 3 - t * ((3. * k ** 2 - 7) / 20)
            poly[:, i] = poly[:, i] / ((poly[:, i] ** 2).sum()) ** .5
        elif i == 3:
            poly[:, i] = t ** 4 - (t ** 2) * ((3 * k ** 2 - 13) /
                                              14.) + 3. * (k ** 2 - 1) * (k ** 2 - 9) / 560
            poly[:, i] = poly[:, i] / ((poly[:, i] ** 2).sum()) ** .5
    printInfo("The polynomial is listed.")
    for i in xrange(degree):
        if k > 6:
            st = ""
            for j in xrange(3):
                st += " {0:>7.4f}".format(poly[j, i])
            st += ' ...'
            for j in xrange(-3, 0):
                st += " {0:>7.4f}".format(poly[j, i])
        else:
            st = ""
            for j in xrange(int(k)):
                st += " {0:>7.4f}".format(poly[j, i])
        printInfo("\t{0:^5s}:{1}".format(label[i], st))
    result = performRegression(poly, y, label=label, **kwargs)
    result['poly'] = poly
    return result
