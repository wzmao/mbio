# -*- coding: utf-8 -*-
"""This module contains some functions for EM analysis.
"""

__author__ = 'Wenzhi Mao'
__all__ = ['genPvalue', 'calcPcutoff', 'showPcutoff', 'transCylinder',
           'showMRCConnection', 'showMRCConnectionEach', 'mrcSegment']


def interpolationball(matrix, index, step, r, **kwargs):
    """Interpolation the value by the radius(ball).
    The Inverse distance weighting is used to weight each value."""

    from numpy import array, arange, floor, ceil

    position = index * step
    w = []
    v = []
    for i in arange(ceil(index[0] - (r / step[0])) // 1, floor(index[0] + (r / step[0])) // 1 + 1):
        for j in arange(ceil(index[1] - (r / step[1])) // 1, floor(index[1] + (r / step[1])) // 1 + 1):
            for k in arange(ceil(index[2] - (r / step[2])) // 1, floor(index[2] + (r / step[2])) // 1 + 1):
                if (((index[0] - i) * step[0]) ** 2 + ((index[1] - j) * step[1]) ** 2 + ((index[2] - k) * step[2]) ** 2) <= r ** 2:
                    w.append(1.0 / ((((index[0] - i) * step[0]) ** 2 + (
                        (index[1] - j) * step[1]) ** 2 + ((index[2] - k) * step[2]) ** 2) ** 2) ** .5)
                    v.append(matrix[i, j, k])
    w = array(w)
    v = array(v)
    w = w / w.sum()
    return (w * v).sum()


def interpolationcube(m, p, way, *kwargs):
    """Interpolation the value by the smallest box.
    The Inverse distance weighting or Trilinear interpolation is
    used to weight each value."""

    from numpy import array

    if way == 'idw':
        tt = array([[[0, 0], [0, 0]], [[0, 0], [0, 0]]], dtype=float)
        tt[0,:,:] += p[0] ** 2
        tt[1,:,:] += (1 - p[0]) ** 2
        tt[:, 0,:] += p[1] ** 2
        tt[:, 1,:] += (1 - p[1]) ** 2
        tt[:,:, 0] += p[2] ** 2
        tt[:,:, 1] += (1 - p[2]) ** 2
        tt = tt ** .5
        tt = 1. / tt
        tt = tt / tt.sum()
    elif way == 'interpolation':
        tt = array([[[1, 1], [1, 1]], [[1, 1], [1, 1]]], dtype=float)
        tt[0,:,:] *= 1 - p[0]
        tt[1,:,:] *= p[0]
        tt[:, 0,:] *= 1 - p[1]
        tt[:, 1,:] *= p[1]
        tt[:,:, 0] *= 1 - p[2]
        tt[:,:, 1] *= p[2]
    if m.shape != (2, 2, 2):
        return -999999999999
    else:
        return (tt * m).sum()


def genPvalue(
    pdb, mrc, sample=None, method=('cube', 'interpolation'), sampleradius=3.0,
        assumenorm=False, **kwargs):
    """This function assign p-value for pdb structure in mrc. sample is used to get the population values.
    `method` must be a tuple or list.
    There are 2 methods now: `cube` and `ball`.
    For the `cube` method, you should provide either ('cube','interpolation') or ('cube','idw').
        `idw` stand for `Inverse distance weighting`, and it is the default option.
    For the `ball` method, you should provide the radius(A) like ('ball',3).
    `sample` should be a `prody.AtomGroup` or a `numpy` n*3 array or `None` to indicate all data to sample.
    `assumenorm` is set to be `False` in default which will get the sample p-value. If it is set
    to `True`, the p-value will be the norm-p-value.
    The p-value will be set to the beta in the pdb.
    """

    from ..IO.output import printError, printInfo
    from .mrc import MRC as MRCclass
    from ..Application.algorithm import binarySearch
    from prody import AtomGroup as pdbclass
    from prody.atomic.selection import Selection as selectionclass
    from numpy import ndarray, zeros_like, array, floor, ceil, rint
    from scipy.stats import norm

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

    if not isinstance(pdb, (pdbclass, selectionclass)):
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
        for i in xrange(int(floor(-sampleradius / step[0])), int(ceil(sampleradius / step[0]) + 1)):
            for j in xrange(int(floor(-sampleradius / step[1])), int(ceil(sampleradius / step[1]) + 1)):
                for k in xrange(int(floor(-sampleradius / step[2])), int(ceil(sampleradius / step[2]) + 1)):
                    if (i * step[0]) ** 2 + (j * step[1]) ** 2 + (k * step[2]) ** 2 <= sampleradius ** 2:
                        ballindex[0].append(i)
                        ballindex[1].append(j)
                        ballindex[2].append(k)
        ballindex = [array(i, dtype=int) for i in ballindex]
        k = array([[len(grid[0])], [len(grid[1])], [len(grid[2])]])
        for i in xrange(len(sample)):
            t = array([ballindex[0] + tempindex[i][0], ballindex[1] +
                       tempindex[i][1], ballindex[2] + tempindex[i][2]])
            t = t[:, (t >= 0).all(0) & (t < k).all()]
            mark[(t[0], t[1], t[2])] = 1
        findset = mrc.data[mark != 0]
    printInfo("Sorting the sample set.")
    findset.sort(kind='quicksort')
    findsetlength = len(findset)
    if assumenorm:
        mu = findset.mean()
        sigma = (
            ((findset - findset.mean()) ** 2).sum() / (findsetlength - 1)) ** .5
    else:
        mu = sigma = 0.

    printInfo("Interpolating the data and assigning p-value.")
    beta = pdb.getBetas()
    coor = pdb.getCoords()
    if method == 'ball':
        index = (coor - gridstart) / step
        for i in xrange(len(coor)):
            beta[i] = interpolationball(mrc.data, index[i], step, r=way)
            if assumenorm:
                beta[i] = norm.cdf(-(beta[i] - mu) / sigma)
            else:
                beta[i] = 1. - \
                    binarySearch(findset, beta[i]) * 1.0 / findsetlength
    elif method == 'cube':
        index = (coor - gridstart) // step
        for i in xrange(len(coor)):
            beta[i] = interpolationcube(mrc.data[index[i][0]:index[i][
                0] + 2, index[i][1]:index[i][1] + 2, index[i][2]:index[i][2] + 2], coor[i] % array(step) / array(step), way)
            if assumenorm:
                beta[i] = norm.cdf(-(beta[i] - mu) / sigma)
            else:
                beta[i] = 1. - \
                    binarySearch(findset, beta[i]) * 1.0 / findsetlength
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


def calcPcutoff(data, scale=5.0, **kwargs):
    """This is a function to calculate the cutoff for high p-values.

        `data` could be a `prody.AtomGroup` with p-values in the Beta, the
    backbone average is calculated to perform analysis. It could also be raw
    number array.

        A linear regression is performed by the first half data and the sigma
    is calculated.
        Cutoff is set to be the first one accepted by `scale`*sigma in the
    tail of data.

        We suggest the cutoff is set for each chain. You need select atom and
    then use this function."""

    from ..IO.output import printError
    from prody import AtomGroup as pdbclass
    from prody.atomic.selection import Selection as selectionclass
    from numpy import ndarray, array, arange

    if isinstance(data, (pdbclass, selectionclass)):
        data = [i for i in data.getHierView().iterResidues()]
        # data.sort(cmp=lambda x,y:chainsort(x.getChid(),y.getChid()) if
        # x.getChid()!=y.getChid() else cmp(x.getResnum(),y.getResnum()))
        data = array([i.select('backbone').getBetas().mean() for i in data])
        data.sort()
    elif isinstance(data, ndarray):
        data.sort()
    else:
        printError(
            "The data format is not supported.(`prody.AtomGroup` or `numpy.ndarray`)")
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
    rmse = (
        ((firsthalfindex * beta1 + beta0 - firsthalfdata) ** 2).sum() / (n - 1)) ** .5
    # Test the second half and get cutoff
    tvalue = abs(index * beta1 + beta0 - data)
    tbigset = (tvalue <= scale * rmse).nonzero()[0]
    return data[max(tbigset)]


def showPcutoff(data, plot, scale=5.0, color=None, detail=False, **kwargs):
    """This is a function to plot the p-value cutoff.

        `data` must be a `prody.AtomGroup` with p-values in the Beta, the
    backbone average is calculated to perform analysis.
        `color` could set to draw in specific color.
        `detail` could be used to plot more detail information.
            1 to plot the error bar. Provide color list.
            2 to also plot sidechain information. Provide 2 colors.

        A linear regression is performed by the first half data and the sigma
    is calculated.
        Cutoff is set to be the first one accepted by `scale`*sigma in the
    tail of data.

        We suggest the cutoff is set for each chain. You need select atom and
    then use this function."""

    from ..IO.output import printError, printInfo
    from prody import AtomGroup as pdbclass
    from prody.atomic.selection import Selection as selectionclass
    from matplotlib.axes import Axes
    from numpy import ndarray, array, arange

    if isinstance(data, (pdbclass, selectionclass)):
        data = [i for i in data.getHierView().iterResidues()]
        data.sort(cmp=lambda x, y: chainsort(x.getChid(), y.getChid()) if x.getChid(
        ) != y.getChid() else cmp(x.getResnum(), y.getResnum()))
        labelindex = array([i.getResnum() for i in data])
    else:
        printError("The data format is not supported.(`prody.AtomGroup`)")
        return None

    data1 = array([i.select('backbone').getBetas().mean() for i in data])
    data1.sort()
    index = arange(len(data1))
    firsthalfdata = array(data1[:len(data1) // 2])
    firsthalfindex = arange(len(firsthalfdata))
    n = len(firsthalfdata)
    # Regression the line
    beta1 = ((firsthalfindex * firsthalfdata).sum() - n * firsthalfindex.mean() * firsthalfdata.mean()) / \
        ((firsthalfindex * firsthalfindex).sum() - n *
         firsthalfindex.mean() * firsthalfindex.mean())
    beta0 = firsthalfdata.mean() - beta1 * firsthalfindex.mean()
    # Determine the RMSE
    rmse = (
        ((firsthalfindex * beta1 + beta0 - firsthalfdata) ** 2).sum() / (n - 1)) ** .5
    # Test the second half and get cutoff
    tvalue = abs(index * beta1 + beta0 - data1)
    tbigset = (tvalue <= scale * rmse).nonzero()[0]
    cutoff = data1[max(tbigset)]

    if isinstance(plot, Axes):
        if detail <= 0:
            if type(color) != type(None):
                plot.plot(labelindex, array(
                    [i.select('backbone').getBetas().mean() for i in data]), '-', c=color, zorder=10)
            else:
                plot.plot(labelindex, array(
                    [i.select('backbone').getBetas().mean() for i in data]), '-', zorder=10)
            plot.plot(
                list(plot.get_xlim()), [cutoff, cutoff], '--', c='grey', alpha=0.8, zorder=5)
            dd = array([i.select('backbone').getBetas().mean() for i in data])
            for i in xrange(len(dd) - 2):
                if (dd[i:i + 3] > cutoff).all():
                    plot.plot(
                        labelindex[i:i + 3], dd[i:i + 3], '.-', c='red', zorder=11)
                elif dd[i] > cutoff:
                    plot.plot(
                        labelindex[i:i + 1], dd[i:i + 1], '.-', c='red', zorder=11)
            x = plot.get_xlim()
            y = plot.get_ylim()
            if detail != -1:
                plot.text(
                    x[1] - 0.05 * (x[1] - x[0]), y[1] - 0.05 * (y[1] - y[0]),
                    "cutoff=%.3f" % (cutoff), va='top', multialignment='left', ha='right')
        elif detail > 0:
            if detail == 2:
                dd = array([i.select('not backbone').getBetas().mean()
                            for i in data])
                yerr3 = dd - \
                    array([i.select('not backbone').getBetas().min()
                           for i in data])
                yerr4 = array(
                    [i.select('not backbone').getBetas().max() for i in data]) - dd
                if type(color) != type(None):
                    plot.plot(
                        labelindex, dd, '-', c=color[1], zorder=10, alpha=0.5)
                    plot.errorbar(labelindex, dd, yerr=[
                                  yerr3, yerr4], capsize=0, elinewidth=0.2, c=color[1], zorder=1)
                else:
                    plot.plot(labelindex, dd, '-', zorder=10, alpha=0.5)
                    plot.errorbar(
                        labelindex, dd, yerr=[yerr3, yerr4], capsize=0, elinewidth=0.1, zorder=1)
            if type(color) != type(None):
                plot.plot(labelindex, array(
                    [i.select('backbone').getBetas().mean() for i in data]), '-', c=color[0], zorder=10)
            else:
                plot.plot(labelindex, array(
                    [i.select('backbone').getBetas().mean() for i in data]), '-', zorder=10)
            plot.plot(
                list(plot.get_xlim()), [cutoff, cutoff], '--', c='grey', alpha=0.8, zorder=5)
            dd = array([i.select('backbone').getBetas().mean() for i in data])
            for i in xrange(len(dd) - 2):
                if (dd[i:i + 3] > cutoff).all():
                    plot.plot(
                        labelindex[i:i + 3], dd[i:i + 3], '.-', c='red', zorder=11)
                elif dd[i] > cutoff:
                    plot.plot(
                        labelindex[i:i + 1], dd[i:i + 1], '.-', c='red', zorder=11)
            yerr1 = dd - \
                array([i.select('backbone').getBetas().min() for i in data])
            yerr2 = array([i.select('backbone').getBetas().max()
                           for i in data]) - dd
            if type(color) != type(None):
                plot.errorbar(labelindex, dd, yerr=[
                              yerr1, yerr2], capsize=0, elinewidth=0.2, c=color[0], zorder=2)
            else:
                plot.errorbar(
                    labelindex, dd, yerr=[yerr1, yerr2], capsize=0, elinewidth=0.2, zorder=2)
            x = plot.get_xlim()
            y = plot.get_ylim()
            plot.text(x[1] - 0.05 * (x[1] - x[0]), y[1] - 0.05 * (y[1] - y[0]),
                      "cutoff=%.3f" % (cutoff), va='top', multialignment='left', ha='right')
    else:
        try:
            if detail <= 0:
                if type(color) != type(None):
                    plot[0].plot(labelindex, array(
                        [i.select('backbone').getBetas().mean() for i in data]), '-', c=color, zorder=10)
                    plot[1].plot(data1, '-', c=color, zorder=10)
                else:
                    plot[0].plot(labelindex, array(
                        [i.select('backbone').getBetas().mean() for i in data]), '-', zorder=10)
                    plot[1].plot(data1, '-', zorder=10)
                plot[0].plot(
                    list(plot[0].get_xlim()), [cutoff, cutoff], '--', c='grey', alpha=0.8, zorder=5)
                dd = array([i.select('backbone').getBetas().mean()
                            for i in data])
                for i in xrange(len(dd) - 2):
                    if (dd[i:i + 3] > cutoff).all():
                        plot[0].plot(
                            labelindex[i:i + 3], dd[i:i + 3], '.-', c='red', zorder=11)
                    elif dd[i] > cutoff:
                        plot[0].plot(
                            labelindex[i:i + 1], dd[i:i + 1], '.-', c='red', zorder=11)
                plot[1].plot(
                    plot[1].get_xlim(), [cutoff, cutoff], '--', c='grey', alpha=0.8, zorder=5)
                plot[1].set_ylim(plot[0].get_ylim())
                x = plot[1].get_xlim()
                y = plot[1].get_ylim()
                if detail != -1:
                    plot[1].text(
                        x[1] - 0.05 * (x[1] - x[0]), y[
                            1] - 0.05 * (y[1] - y[0]),
                        "cutoff=%.3f" % (cutoff), va='top', multialignment='left', ha='right')
            elif detail > 0:
                if detail == 2:
                    dd = array(
                        [i.select('not backbone').getBetas().mean() for i in data])
                    yerr3 = dd - \
                        array([i.select('not backbone').getBetas().min()
                               for i in data])
                    yerr4 = array(
                        [i.select('not backbone').getBetas().max() for i in data]) - dd
                    if type(color) != type(None):
                        plot[0].plot(
                            labelindex, dd, '-', c=color[1], zorder=10, alpha=0.5)
                        plot[0].errorbar(labelindex, dd, yerr=[
                                         yerr3, yerr4], capsize=0, elinewidth=0.2, c=color[1], zorder=1)
                    else:
                        plot[0].plot(labelindex, dd, '-', zorder=10, alpha=0.5)
                        plot[0].errorbar(
                            labelindex, dd, yerr=[yerr3, yerr4], capsize=0, elinewidth=0.1, zorder=1)
                if type(color) != type(None):
                    plot[0].plot(labelindex, array(
                        [i.select('backbone').getBetas().mean() for i in data]), '-', c=color[0], zorder=10)
                    plot[1].plot(data1, '-', c=color[0], zorder=10)
                else:
                    plot[0].plot(labelindex, array(
                        [i.select('backbone').getBetas().mean() for i in data]), '-', zorder=10)
                    plot[1].plot(data1, '-', zorder=10)
                plot[0].plot(
                    list(plot[0].get_xlim()), [cutoff, cutoff], '--', c='grey', alpha=0.8, zorder=5)
                dd = array([i.select('backbone').getBetas().mean()
                            for i in data])
                for i in xrange(len(dd) - 2):
                    if (dd[i:i + 3] > cutoff).all():
                        plot[0].plot(
                            labelindex[i:i + 3], dd[i:i + 3], '.-', c='red', zorder=11)
                    elif dd[i] > cutoff:
                        plot[0].plot(
                            labelindex[i:i + 1], dd[i:i + 1], '.-', c='red', zorder=11)
                yerr1 = dd - \
                    array([i.select('backbone').getBetas().min()
                           for i in data])
                yerr2 = array([i.select('backbone').getBetas().max()
                               for i in data]) - dd
                if type(color) != type(None):
                    plot[0].errorbar(labelindex, dd, yerr=[
                                     yerr1, yerr2], capsize=0, elinewidth=0.2, c=color[0], zorder=2)
                else:
                    plot[0].errorbar(
                        labelindex, dd, yerr=[yerr1, yerr2], capsize=0, elinewidth=0.2, zorder=2)
                plot[1].plot(
                    plot[1].get_xlim(), [cutoff, cutoff], '--', c='grey', alpha=0.8, zorder=5)
                plot[1].set_ylim(plot[0].get_ylim())
                x = plot[1].get_xlim()
                y = plot[1].get_ylim()
                plot[1].text(
                    x[1] - 0.05 * (x[1] - x[0]), y[1] - 0.05 * (y[1] - y[0]),
                    "cutoff=%.3f" % (cutoff), va='top', multialignment='left', ha='right')
        except:
            printError(
                "The plot type wrong. Must be 1 or 2 `matplotlib.axes.Axes`.")
            return None


def genPvalueSample(mrc, sample=None, sampleradius=3.0, **kwargs):
    """Given the `mrc` and a sample structure, return the sample set around the sample
    structure with radius `sampleradius`.
    """

    from ..IO.output import printError, printInfo
    from .mrc import MRC as MRCclass
    from prody import AtomGroup as pdbclass
    from numpy import ndarray, zeros_like, array, floor, ceil, rint

    if not isinstance(mrc, MRCclass):
        printError("Only mbio.MRC class supported for `mrc`.")
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
        for i in xrange(int(floor(-sampleradius / step[0])), int(ceil(sampleradius / step[0]) + 1)):
            for j in xrange(int(floor(-sampleradius / step[1])), int(ceil(sampleradius / step[1]) + 1)):
                for k in xrange(int(floor(-sampleradius / step[2])), int(ceil(sampleradius / step[2]) + 1)):
                    if (i * step[0]) ** 2 + (j * step[1]) ** 2 + (k * step[2]) ** 2 <= sampleradius ** 2:
                        ballindex[0].append(i)
                        ballindex[1].append(j)
                        ballindex[2].append(k)
        ballindex = [array(i, dtype=int) for i in ballindex]
        k = array([[len(grid[0])], [len(grid[1])], [len(grid[2])]])
        for i in xrange(len(sample)):
            t = array([ballindex[0] + tempindex[i][0], ballindex[1] +
                       tempindex[i][1], ballindex[2] + tempindex[i][2]])
            t = t[:, (t >= 0).all(0) & (t < k).all()]
            mark[(t[0], t[1], t[2])] = 1
        findset = mrc.data[mark != 0]
    printInfo("Sorting the sample set.")
    findset.sort(kind='quicksort')
    return findset


def transCylinder(pdb, **kwargs):
    """Transfer the PDB fragment to Cylinder.
    Given the PDB, extract the CA C N and generate the center, 3 directions of the cube
    and the length."""

    from ..IO.output import printInfo
    from numpy.linalg import svd
    from numpy import cross, array

    p1 = pdb.select("name CA C N").copy()
    if p1 is None:
        printError("The pdb has no CA C or N in the backbone.")
        return None

    if p1.numAtoms() < 15:
        printError("The atom number({0}) is not enough to perform calculation.".format(
            p1.numAtoms()))
        printError("The result is not reliable.")

    data = p1.getCoords()
    datamean = data.mean(axis=0)
    uu, dd, vv = svd(data - datamean)
    cent = datamean
    dirc1 = vv[0]
    if abs((dirc1 ** 2).sum() ** .5 - 1) > 1e-10:
        raise ValueError(
            "length of dirc is not 1, is {0}".format((dirc1 ** 2).sum() ** .5))
    if (data[-1] - data[0]).dot(dirc1) < 0:
        dirc1 = -dirc1
    rank = (data - cent).dot(dirc1)
    rankrange = [rank.min(), rank.max()]
    dirc2 = cent - cent.dot(dirc1) * dirc1
    dirc2 = dirc2 / ((dirc2 ** 2).sum() ** .5)
    dirc3 = cross(dirc1, dirc2, axis=0)
    dirc3 = dirc3 / ((dirc3 ** 2).sum() ** .5)
    dirc = array((dirc1, dirc2, dirc3))
    return cent, dirc, rankrange

def showMRCConnection(mrc, cutoff=2, **kwargs):
    """Plot 3D plot of connected parts in different color for MRC."""

    from matplotlib import pyplot as plt
    from matplotlib import use as matplotlibuse
    import mpl_toolkits.mplot3d.axes3d as p3
    from mpl_toolkits.mplot3d import Axes3D
    from ..Application.setting import getMatplotlibDisplay
    from ..Application.plotting import setAxesEqual
    from numpy import array

    try:
        if not getMatplotlibDisplay():
            matplotlibuse('Agg')
    except:
        pass

    fig = plt.figure(figsize=(6, 6), facecolor='white')
    ax = p3.Axes3D(fig, aspect=1)
    ax.w_xaxis.set_pane_color((0, 0, 0))
    ax.w_yaxis.set_pane_color((0, 0, 0))
    ax.w_zaxis.set_pane_color((0, 0, 0))
    ax.w_xaxis.line.set_lw(0)
    ax.w_yaxis.line.set_lw(0)
    ax.w_zaxis.line.set_lw(0)
    classes = {}
    step = mrc.getGridSteps()
    cutoff = cutoff**2
    for i, j, k in zip(*mrc.data.nonzero()):
        if mrc.data[i, j, k] == 0:
            continue
        if mrc.data[i, j, k] not in classes.keys():
            classes[mrc.data[i, j, k]] = [[i, j, k]]
        else:
            classes[mrc.data[i, j, k]].append([i, j, k])
    for ty, i in zip(classes.keys(), xrange(len(classes))):
        color = plt.cm.gist_ncar(i*1./len(classes)*.9)
        pos = array(classes[ty])
        ax.scatter(pos[:, 0]*step[0]+mrc.origin[0],
                   pos[:, 1]*step[1]+mrc.origin[1],
                   pos[:, 2]*step[2]+mrc.origin[2], lw=0, c=color, zorder=10)
        if cutoff>0:
            for j in xrange(len(pos)):
                for k in xrange(j):
                    if (((pos[j]-pos[k])*step)**2).sum() <= cutoff:
                        ax.plot(pos[[j, k], 0]*step[0]+mrc.origin[0],
                                pos[[j, k], 1]*step[1]+mrc.origin[1],
                                pos[[j, k], 2]*step[2]+mrc.origin[2], lw=3, c=color, zorder=10)
        del pos
    ax.set_xlabel('X', fontsize=15)
    ax.set_ylabel('Y', fontsize=15)
    ax.set_zlabel('Z', fontsize=15)
    setAxesEqual(ax)
    del classes
    del ax
    plt.ion()
    try:
        if getMatplotlibDisplay():
            plt.show()
    except:
        pass
    return fig

def showMRCConnectionEach(mrc, cutoff=2, path=None, **kwargs):
    """Plot 3D plot of connected parts in different color for MRC."""

    from matplotlib import pyplot as plt
    from matplotlib import use as matplotlibuse
    import mpl_toolkits.mplot3d.axes3d as p3
    from mpl_toolkits.mplot3d import Axes3D
    from ..Application.setting import getMatplotlibDisplay
    from ..Application.plotting import setAxesEqual
    from numpy import array
    import os

    try:
        if not getMatplotlibDisplay():
            matplotlibuse('Agg')
    except:
        pass

    if path is None:
        path=os.getcwd()
    fig = plt.figure(figsize=(6, 6), facecolor='white')
    ax = p3.Axes3D(fig, aspect=1)
    ax.w_xaxis.set_pane_color((0, 0, 0))
    ax.w_yaxis.set_pane_color((0, 0, 0))
    ax.w_zaxis.set_pane_color((0, 0, 0))
    ax.w_xaxis.line.set_lw(0)
    ax.w_yaxis.line.set_lw(0)
    ax.w_zaxis.line.set_lw(0)
    classes = {}
    step = mrc.getGridSteps()
    grid = mrc.getGridCoords()
    cutoff = cutoff**2
    for i, j, k in zip(*mrc.data.nonzero()):
        if mrc.data[i, j, k] == 0:
            continue
        if mrc.data[i, j, k] not in classes.keys():
            classes[mrc.data[i, j, k]] = [[i, j, k]]
        else:
            classes[mrc.data[i, j, k]].append([i, j, k])
    sca=ax.scatter([60,240],[60,240],[60,240],lw=0,zorder=10)
    plt.ion()
    ax.set_xlabel('X', fontsize=15)
    ax.set_ylabel('Y', fontsize=15)
    ax.set_zlabel('Z', fontsize=15)
    setAxesEqual(ax)
    for ty, i in zip(classes.keys(), xrange(len(classes))):
        color = plt.cm.gist_ncar(i*1./len(classes)*.9)
        pos = array(classes[ty])
        sca._offsets3d=pos[:, 0]*step[0]+mrc.origin[0],pos[:, 1]*step[1]+mrc.origin[1],pos[:, 2]*step[2]+mrc.origin[2]
        sca._facecolor3d=color
        del pos
        plt.savefig(os.path.join(path,str(i)+'.png'))
    del classes
    del ax
    return fig

def testfit(pos,step):

    from numpy import array, diag
    from numpy.linalg import svd
    
    data=array(pos)*step
    datamean = data.mean(axis=0)
    uu, dd, vv = svd(data - datamean,full_matrices=False)
    d=dd**2
    dd[0]=0
    if (((uu.dot(diag(dd)).dot(vv))**2).sum(1)**.5<3).all():
        return 1
    else:
        return 2

    # Old 3
    # if (((uu.dot(diag(dd)).dot(vv))**2).sum(1)**.5<6).sum()*1./data.shape[0]>.9:
    #     return 1
    # else:
    #     return 2

    # Old 2
    # if d[0]/d.sum() <.8:
    #     return 2
    # else:
    #     return 1

    # Old 1
    # if len(pos)>30:
    #     return 2
    # else:
    #     return 1

def mrcSegment(mrc, percentage=0.001, cutoff=3, **kwargs):
    """Segment the MRC with the top `percentage` points.
    Only two points closer than the cutoff will be taken as connected."""

    from numpy import floor, ceil, argsort, zeros, zeros_like, array, unique
    from ..IO.output import printUpdateInfo, finishUpdate
    from .mrc import MRC

    maxnum = int(percentage*mrc.data.size)
    args = argsort(mrc.data.ravel())[:mrc.data.size-maxnum-1:-1]
    pos = zeros((maxnum, 3), dtype=int)
    pos[:, 0] = args // (mrc.data.shape[1]*mrc.data.shape[2])
    pos[:, 1] = args % (mrc.data.shape[1]*mrc.data.shape[2]) // mrc.data.shape[2]
    pos[:, 2] = args % (mrc.data.shape[2])
    data = mrc.data.ravel()[args]

    save=zeros_like(mrc.data,dtype=int)
    save[[pos[:,0],pos[:,1],pos[:,2]]]=-1
    save1=zeros_like(mrc.data,dtype=float)

    origin = mrc.origin
    grid=mrc.getGridSteps()
    step = cutoff/grid
    ranges=[xrange(int(floor(-step[i])),int(ceil(step[i])+int(step[i].is_integer()))) for i in xrange(3)]
    cutoff2=cutoff**2

    classnum=0
    save1count=0
    classmp={}
    classmpreverse={}
    classcount={}
    classpos={}

    for posnum in xrange(maxnum):
        if posnum%1000==0:
            printUpdateInfo("Building {:10d}/{:10d}".format(posnum,maxnum))
        temp=pos[posnum]
        closeset=[]
        closetype=[]
        closenumber=0
        for i in ranges[0]:
            for j in ranges[1]:
                for k in ranges[2]:
                    if save[temp[0]+i,temp[1]+j,temp[2]+k]>0 and (i*grid[0])**2+(j*grid[1])**2+(k*grid[2])**2<=cutoff2:
                        closeset.append([temp+array([i,j,k])])
                        closetype.append(classmp[save[temp[0]+i,temp[1]+j,temp[2]+k]])
                        closenumber+=1
        if closenumber==0:
            classnum+=1
            save[temp[0],temp[1],temp[2]]=classnum
            classcount[classnum]=[1,0]
            classmp[classnum]=classnum
            classmpreverse[classnum]=[classnum]
            classpos[classnum]=[pos[posnum]]
        elif len(unique(closetype))==1:
            typeclass=closetype[0]
            save[temp[0],temp[1],temp[2]]=typeclass
            orilen=classcount[typeclass][0]
            classcount[typeclass][0]+=1
            classpos[typeclass].append(pos[posnum])
            if classcount[typeclass][1]==0:
                if classcount[typeclass][0]>=10:
                    classcount[typeclass][1]=testfit(classpos[typeclass],grid)
                    if classcount[typeclass][1]==2:
                        save1count+=1
                        tempposlist = classpos[typeclass]
                        for i in xrange(orilen):
                            save1[tempposlist[i][0],tempposlist[i][1],tempposlist[i][2]]=save1count
                        del tempposlist
            elif classcount[typeclass][1]==1:
                classcount[typeclass][1]=testfit(classpos[typeclass],grid)
                if classcount[typeclass][1]==2:
                    save1count+=1
                    tempposlist = classpos[typeclass]
                    for i in xrange(orilen):
                        save1[tempposlist[i][0],tempposlist[i][1],tempposlist[i][2]]=save1count
                    del tempposlist
            else:
                pass
            del typeclass
        else:
            closetypesort=unique(closetype)
            typeclass=closetypesort[0]
            save[temp[0],temp[1],temp[2]]=typeclass
            orilen=classcount[typeclass][0]
            classcount[typeclass][0]+=1
            classpos[typeclass].append(pos[posnum])
            hasnocylinder=False
            for i in closetypesort[1:]:
                if classcount[i][1]==2:
                    hasnocylinder=True
                classcount[typeclass][0]+=classcount[i][0]
                classpos[typeclass]+=classpos[i]
                classmp[i]=typeclass
                for j in classmpreverse[i]:
                    classmp[j] = typeclass
                classmpreverse[typeclass]+=classmpreverse[i]
                classmpreverse.pop(i)
            if classcount[typeclass][1]==0:
                if classcount[typeclass][0]>=10:
                    classcount[typeclass][1]=testfit(classpos[typeclass],grid) if not hasnocylinder else 2
                    if classcount[typeclass][1]==2:
                        for i in closetypesort[1:]:
                            if classcount[i][1]==1:
                                save1count+=1
                                tempposlist = classpos[i]
                                for i in xrange(len(classpos[i])):
                                    save1[tempposlist[i][0],tempposlist[i][1],tempposlist[i][2]]=save1count
                                del tempposlist
            elif classcount[typeclass][1]==1:
                classcount[typeclass][1]=testfit(classpos[typeclass],grid) if not hasnocylinder else 2
                if classcount[typeclass][1]==2:
                    for i in closetypesort[1:]:
                        if classcount[i][1]==1:
                            save1count+=1
                            tempposlist = classpos[i]
                            for i in xrange(len(classpos[i])):
                                save1[tempposlist[i][0],tempposlist[i][1],tempposlist[i][2]]=save1count
                            del tempposlist
                    save1count+=1
                    tempposlist = classpos[typeclass]
                    for i in xrange(orilen):
                        save1[tempposlist[i][0],tempposlist[i][1],tempposlist[i][2]]=save1count
                    del tempposlist
            else:
                pass
            for i in closetypesort[1:]:
                classcount.pop(i)
                classpos.pop(i)
            del typeclass,closetypesort
        del temp,closeset,closetype,closenumber

    for i in classcount:
        if classcount[i][1]==1:
            save1count+=1
            tempposlist = classpos[i]
            for i in xrange(len(tempposlist)):
                save1[tempposlist[i][0],tempposlist[i][1],tempposlist[i][2]]=save1count

    del classnum,save1count,classmp,classmpreverse,classcount,classpos

    finishUpdate()
    mrc1=MRC()
    for i in mrc.header.__dict__:
        setattr(mrc1.header,i,getattr(mrc.header,i))

    mrc1.data=save1
    mrc1.update()

    return mrc1
