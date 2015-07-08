# -*- coding: utf-8 -*-
"""This module contains some functions for EM analysis.
"""

__author__ = 'Wenzhi Mao'
__all__ = ['genPvalue', 'calcPcutoff']


def interpolationball(matrix, index, step, r, **kwarg):
    """Interpolation the value by the radius(ball).
    The Inverse distance weighting is used to weight each value."""

    from numpy import array, arange, floor, ceil

    position = index * step
    w = []
    v = []
    for i in arange(ceil(index[0] - (r / step[0])) // 1, floor(index[0] + (r / step[0])) // 1 + 1):
        for j in arange(ceil(index[1] - (r / step[1])) // 1, floor(index[1] + (r / step[1])) // 1 + 1):
            for k in arange(ceil(index[2] - (r / step[2])) // 1, floor(index[2] + (r / step[2])) // 1 + 1):
                if (((index[0] - i) * step[0])**2 + ((index[1] - j) * step[1])**2 + ((index[2] - k) * step[2])**2) <= r**2:
                    w.append(1.0 / ((((index[0] - i) * step[0])**2 + (
                        (index[1] - j) * step[1])**2 + ((index[2] - k) * step[2])**2)**2)**.5)
                    v.append(matrix[i, j, k])
    w = array(w)
    v = array(v)
    w = w / w.sum()
    return (w * v).sum()


def interpolationcube(m, p, way, *kwarg):
    """Interpolation the value by the smallest box.
    The Inverse distance weighting or Trilinear interpolation is
    used to weight each value."""

    from numpy import array

    if way == 'idw':
        tt = array([[[0, 0], [0, 0]], [[0, 0], [0, 0]]], dtype=float)
        tt[0, :, :] += p[0]**2
        tt[1, :, :] += (1 - p[0])**2
        tt[:, 0, :] += p[1]**2
        tt[:, 1, :] += (1 - p[1])**2
        tt[:, :, 0] += p[2]**2
        tt[:, :, 1] += (1 - p[2])**2
        tt = tt**.5
        tt = 1. / tt
        tt = tt / tt.sum()
    elif way == 'interpolation':
        tt = array([[[1, 1], [1, 1]], [[1, 1], [1, 1]]], dtype=float)
        tt[0, :, :] *= 1 - p[0]
        tt[1, :, :] *= p[0]
        tt[:, 0, :] *= 1 - p[1]
        tt[:, 1, :] *= p[1]
        tt[:, :, 0] *= 1 - p[2]
        tt[:, :, 1] *= p[2]
    return (tt * m).sum()

def genPvalue(pdb, mrc, sample=None, method=('cube', 'interpolation'), sampleradius=3.0,
              **kwarg):
    """`method` must be a tuple or list.
    There are 2 methods now: `cube` and `ball`.
    For the `cube` method, you should provide either ('cube','interpolation') or ('cube','idw').
        `idw` stand for `Inverse distance weighting`, and it is the default option.
    For the `ball` method, you should provide the radius(A) like ('ball',3).
    `sample` should be a `prody.AtomGroup` or a `numpy` n*3 array or `None` to indicate all data to sample.
    The p-value will be set to the beta in the pdb.
    """

    from ..IO.output import printError, printInfo
    from ..IO.mrc import MRC as MRCclass
    from ..Application.algorithm import binarySearch
    from prody import AtomGroup as pdbclass
    from numpy import ndarray, zeros_like, array, floor, ceil, rint

    if isinstance(method, (tuple, list)):
        if len(method) == 0:
            printError("The method is not valid.")
            return None
        if method[0].lower() == 'cube':
            if len(method) == 1:
                way = 'idw'
                method = 'cube'
            elif method[1].lower() in ['interpolation', 'idw']:
                way = method[1].lower()
                method = 'cube'
            else:
                printError("The method[1] is not valid.")
                printError(
                    "Only 'idw' or 'interpolation' supported for 'cube'.")
                return None
        elif method[0].lower() == 'ball':
            if len(method) < 2:
                printError("The radius must provide for the `ball` method.")
            try:
                way = float(eval(str(method[1])))
            except:
                printError(
                    "Only numbers are support as the second option for `ball`.")
                return None
            if way <= 0:
                printError("Radius must be positive.")
                return None
            else:
                method = 'ball'
    elif isinstance(method, (str)):
        if method.lower() == 'cube':
            method = 'cube'
            way = 'idw'
        else:
            printError("Only `cube` support no option format.")
            return None
    else:
        printError("The method must be tuple or list")
        return None

    if not isinstance(mrc, MRCclass):
        printError("Only mbio.MRC class supported for `mrc`.")
        return None

    if not isinstance(pdb, pdbclass):
        printError("Only prody.AtomGroup class supported for `pdb`.")
        return None

    if type(sample) == type(None):
        sample = None
    elif isinstance(sample, ndarray):
        if not (sample.shape == (3,) or (len(sample.shape) == 2 and sample.shape[1] == 3)):
            printError("The sample coordinates must has 3 columns.")
            return None
        if sample.shape == (3,):
            sample = array([sample])
    elif isinstance(sample, pdbclass):
        sample = sample.getCoords()

    printInfo("Getting the sample set.")
    mark = zeros_like(mrc.data)
    grid = array(mrc.getGridCoords())
    gridstart = array([grid[0, 0], grid[1, 0], grid[2, 0]])
    step = mrc.getGridSteps()
    if type(sample) == type(None):
        findset = mrc.data.flatten()
    else:
        tempindex = array(
            rint(array(((sample - grid[:, 0]) / step), dtype=float)), dtype=int)
        ballindex = ([], [], [])
        for i in range(int(floor(-sampleradius / step[0])), int(ceil(sampleradius / step[0]) + 1)):
            for j in range(int(floor(-sampleradius / step[1])), int(ceil(sampleradius / step[1]) + 1)):
                for k in range(int(floor(-sampleradius / step[2])), int(ceil(sampleradius / step[2]) + 1)):
                    if (i * step[0])**2 + (j * step[1])**2 + (k * step[2])**2 <= sampleradius**2:
                        ballindex[0].append(i)
                        ballindex[1].append(j)
                        ballindex[2].append(k)
        ballindex = [array(i, dtype=int) for i in ballindex]
        k = array([[len(grid[0])], [len(grid[1])], [len(grid[2])]])
        for i in range(len(sample)):
            t = array([ballindex[0] + tempindex[i][0], ballindex[1] +
                       tempindex[i][1], ballindex[2] + tempindex[i][2]])
            t = t[:, (t >= 0).all(0) & (t < k).all()]
            mark[(t[0], t[1], t[2])] = 1
        findset = mrc.data[mark != 0]
    printInfo("Sorting the sample set.")
    findset.sort(kind='quicksort')
    findsetlength = len(findset)

    printInfo("Interpolating the data and assigning p-value.")
    beta = pdb.getBetas()
    coor = pdb.getCoords()
    if method == 'ball':
        index = (coor - gridstart) / step
        for i in range(len(coor)):
            beta[i] = interpolationball(mrc.data, index[i], step, r=way)
            beta[i] = 1. - binarySearch(findset, beta[i]) * 1.0 / findsetlength
    elif method == 'cube':
        index = (coor - gridstart) // step
        for i in range(len(coor)):
            beta[i] = interpolationcube(mrc.data[index[i][0]:index[i][
                0] + 2, index[i][1]:index[i][1] + 2, index[i][2]:index[i][2] + 2], coor[i] % array(step) / array(step), way)
            beta[i] = 1. - binarySearch(findset, beta[i]) * 1.0 / findsetlength
    pdb.setBetas(beta)
    return pdb

def chainsort(x, y):
    """Chain id sort function. A-Z then number."""
    if x == y:
        return cmp(0, 0)
    elif x.isdigit() == y.isdigit() == True:
        return cmp(float(x), float(y))
    elif x.isdigit() == y.isdigit() == False:
        return cmp(x, y)
    elif x.isdigit():
        return cmp(2, 1)
    else:
        return cmp(1, 2)

def calcPcutoff(data, scale=5.0, **kwarg):
    """This is a function to calculate the cutoff for high p-values.

        `data` would be a `prody.AtomGroup` with p-values in the Beta, the
    backbone average is calculated to perform analysis. It could also be raw
    number array.

        A linear regression is performed by the first half data and the sigma
    is calculated.
        Cutoff is set to be the first one accepted by `scale`*sigma in the
    tail of data.

        We suggest the cutoff is set for each chain. You need select atom and
    then use this function."""

    from ..IO.output import printError, printInfo
    from prody import AtomGroup as pdbclass
    from numpy import ndarray, array, arange

    if isinstance(data,pdbclass):
        data=[i for i in data.iterResidues()]
        # data.sort(cmp=lambda x,y:chainsort(x.getChid(),y.getChid()) if x.getChid()!=y.getChid() else cmp(x.getResnum(),y.getResnum()))
        data=array([i.select('backbone').getBetas().mean() for i in data])
        data.sort()
    elif isinstance(data,ndarray):
        data.sort()
    else:
        printError("The data format is not supported.(`prody.AtomGroup` or `numpy.ndarray`)")
        return None

    index = arange(len(data))
    firsthalfdata = array(data[:len(data) // 2])
    firsthalfindex = arange(len(firsthalfdata))
    n = len(firsthalfdata)
    # Regression the line
    beta1 = ((firsthalfindex * firsthalfdata).sum() - n * firsthalfindex.mean() * firsthalfdata.mean()) / \
        ((firsthalfindex * firsthalfindex).sum() - n *
         firsthalfindex.mean() * firsthalfindex.mean())
    beta0 = firsthalfdata.mean() - beta1 * firsthalfindex.mean()
    # Determine the RMSE
    rmse = (((firsthalfindex * beta1 + beta0 - firsthalfdata)**2).sum()/(n-1))**.5
    # Test the
    print len(index),len(data)
    tvalue = abs(index * beta1 + beta0 - data)
    tbigset = (tvalue <= scale*rmse).nonzero()[0]
    return data[max(tbigset)]