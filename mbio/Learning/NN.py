# -*- coding: utf-8 -*-
"""This module contains some Neural Networks functions.
"""

__author__ = 'Wenzhi Mao'

__all__ = ["NeuralNet"]


class NeuralNet(object):

    '''This is a neu web class to setup a web and perform the training.
    The content is based on the Artificial Neural Networks.
    The training is performed by the BP algorithm.'''

    def __init__(self, webshape=None, label=None, data=None, **kwargs):
        '''the init function of `neu_net` class. webshape, label, data
        could be assigned by options.'''

        self.inputdata = None
        self.outputdata = None
        self.label = None
        self.trans = None
        self.webshape = None

        self._clockdict = {}

        if webshape is not None:
            self.setWebshape(webshape)
        if label is not None:
            self.setLabel(label)
        if data is not None and len(data) == 2:
            self.setInput(data[0])
            self.setOutput(data[1])

    def _tic(self, key, **kwargs):
        '''Add a start time for some key event.'''
        from time import time
        self._clockdict[key] = time()

    def _toc(self, key, **kwargs):
        '''Remove the key and return the time used for the key event.'''
        from time import time
        tic = self._clockdict.pop(key, None)
        if tic:
            return time() - tic
        else:
            return None

    def __str__(self):
        return "Neural Networks from mbio"

    def _clear_data(self, keep=[], **kwargs):
        '''Delete all data in this class.
        You could use keep to keep the variables you want to keep.'''
        if not ('webshape' in keep or 'shape' in keep):
            self.webshape = None
        if not ('input' in keep or 'inputdata' in keep):
            self.inputdata = None
        if not ('output' in keep or 'outputdata' in keep):
            self.outputdata = None
        if not ('label' in keep):
            self.label = None
        if not ('trans' in keep):
            self.trans = None
        self._clockdict = {}

    def setWebshape(self, webshape=[], **kwargs):
        '''Set the web shape (nodes for each layer) for the network.'''

        from ..IO.output import printInfo, printError
        from numpy import array, int32, zeros, float64, matrix, random, ndarray

        if len(webshape) == 0:
            printError('No web shape data provided. Ignored.')
            return None
        webshape = array(webshape, dtype=int32)
        if isinstance(self.webshape, ndarray) and (self.webshape == webshape).all():
            printInfo('webshape unchanged, the trans also unchanged.', 1)
        else:
            if self.inputdata is not None and isinstance(self.inputdata, ndarray) and self.inputdata.shape[0] != webshape[0]:
                printError("The web shape and the inputdata doesn't fit.")
                self.inputdata = None
                printError("The Input data has been deleted.")
            if self.outputdata is not None and isinstance(self.outputdata, ndarray) and self.outputdata.shape[0] != webshape[-1]:
                printError("The web shape and the outputdata doesn't fit.")
                self.outputdata = None
                printError("The Output data has been deleted.")
            self.webshape = webshape
            self.trans = zeros((webshape.shape[0] - 1), dtype=ndarray)
            for i in range(len(webshape) - 1):
                self.trans[i] = matrix(
                    random.random((webshape[i] + 1, webshape[i + 1])) * 0.1 - 0.05, dtype=float64)
            self.label = ["x{0}".format(i + 1) for i in xrange(webshape[0])]
        return None

    def setLabel(self, label=[], **kwargs):
        '''Set the labels for the inputs.'''

        from ..IO.output import printInfo, printError

        if len(label) == 0:
            printError('No label provided. Ignored.')
            return None
        if not isinstance(label, list):
            label = list(label)
        label = [str(i) for i in label]
        if self.webshape is not None:
            if self.webshape[0] != len(label):
                printError(
                    "The length of the label list doesn't fit the webshape[0].")
            else:
                self.label = label
        else:
            self.label = label

    def setInput(self, inputdata=None, **kwargs):
        '''Set the input data for the network.'''

        from ..IO.output import printInfo, printError
        from numpy import ndarray, matrix, float64, zeros, int32

        if inputdata is None:
            self.printError('No input data provided. Ignored.')
            return None
        if isinstance(inputdata, (ndarray, list)):
            try:
                inputdata = matrix(inputdata, dtype=float64)
            except:
                printError("Input data format mistake.")
                return None
        if isinstance(inputdata, matrix):
            if not inputdata.dtype == float64:
                inputdata = matrix(inputdata, dtype=float64)
        else:
            printError(
                "Please provide a list, numpy.ndarray or a numpy.matrix.")
            return None
        if inputdata.ndim != 2:
            printError("The input must be 2-dimension data.")
            return None
        if self.webshape is not None:
            if inputdata.shape[1] != self.webshape[0]:
                printError(
                    "The input data shape[1] doesn't fit the webshape[0].")
                return None
        if self.outputdata is not None:
            if self.inputdata.shape[0] != outputdata.shape[0]:
                printError("The sizes of Input and Output are different.")
                return None
        self.inputdata = inputdata
        return None

    def setOutput(self, outputdata=None, **kwargs):
        '''Set the output data for the network.'''

        from ..IO.output import printInfo, printError
        from numpy import ndarray, matrix, float64, zeros, int32

        if outputdata is None:
            printError('No output data provided. Ignored.')
            return None
        if isinstance(outputdata, (ndarray, list)):
            try:
                outputdata = matrix(outputdata, dtype=float64)
            except:
                printError("Output data format mistake.")
                return None
        if isinstance(outputdata, matrix):
            if not outputdata.dtype == float64:
                outputdata = matrix(outputdata, dtype=float64)
        else:
            printError(
                "Please provide a list, numpy.ndarray or a numpy.matrix.")
            return None
        if outputdata.ndim != 2:
            printError("The output must be 2-dimension data.")
            return None
        if self.webshape is not None:
            if outputdata.shape[1] != self.webshape[-1]:
                printError(
                    "The output data shape[1] doesn't fit the webshape[-1].")
                return None
        if self.inputdata is not None:
            if self.inputdata.shape[0] != outputdata.shape[0]:
                printError("The sizes of Input and Output are different.")
                return None
        self.outputdata = outputdata
        return None

    def fit(self, times, step='auto', **kwargs):
        '''Train the data.'''

        from .CNN_p import fit_ANN_BP
        from ..IO.output import printInfo, printError

        if self.inputdata.shape[0] != self.outputdata.shape[0]:
            printError("The sizes of Input and Output are different.")
            return None
        if self.inputdata.shape[1] != self.webshape[0]:
            printError("The sizes of Input doesnot fit the webshape.")
            return None
        if self.outputdata.shape[1] != self.webshape[-1]:
            printError("The sizes of Output doesnot fit the webshape.")
            return None
        for i in range(len(self.webshape) - 1):
            if self.trans[i].shape != (self.webshape[i] + 1, self.webshape[i + 1]):
                printError("`trans`[{0}] shape wrong".format(i))
                result = 0
                return None
        result = fit_ANN_BP(shape=self.webshape, input=self.inputdata,
                            output=self.outputdata, trans=self.trans, times=times)
        if isinstance(result, list):
            if result[0] is None:
                printError(result[1])
                return None
        return result

    # def __simulate_point(self, inputd, outputd, step=0.1, **kwargs):
    #     '''Train one data point.'''
    #     outputd1 = np.array(outputd)
    #     layer = len(self.webshape) - 1
    #     temp = inputd
    #     save = [np.array(inputd)]
    #     for i in range(layer):
    #         temp = np.concatenate((temp, np.ones((1, 1))), axis=1)
    #         temp = temp.dot(self.trans[i])
    #         temp = 1. / (1. + np.exp(-temp))
    #         save = save + [np.array(temp)]
    #     wucha = [save[-1] * (1 - save[-1]) * (outputd1 - save[-1])]
    #     for i in reversed(range(layer)):
    #         wucha = [
    #             save[i] * (1 - save[i]) * np.array(wucha[0].dot(self.trans[i][:-1].T))] + wucha
    #     for i in range(layer):
    #         self.trans[i] = self.trans[
    #             i] + (np.concatenate((save[i], np.ones((1, 1))), axis=1).T.dot(wucha[i + 1])) * step
    # return ((outputd1 - save[-1])*(outputd1 - save[-1])).sum()

    def predict(self, inputtest, **kwargs):

        from ..IO.output import printInfo, printError
        from numpy import ndarray, array

        if not isinstance(inputtest, ndarray):
            inputtest = array(inputtest)

        # for i in range(inputtest.shape[0]):
        #     layer = len(self.webshape) - 1
        #     temp = inputtest[i]
        #     for j in range(layer):
        #         temp = np.concatenate((temp, np.ones((1, 1))), axis=1)
        #         temp = temp.dot(self.trans[j])
        #     mark += [f(temp) == f(outputtest[i])]
        # return mark

    def _check(self, **kwargs):
        '''Check all data compatible to each other.
        Check only the data for calculation.'''

        from ..IO.output import printInfo, printError
        from numpy import ndarray, int32, float64

        result = 1
        # Type Check
        # webshape
        if not isinstance(self.webshape, ndarray):
            printInfo("`webshape` type")
            result = 0
        else:
            if self.webshape.dtype != int32:
                printInfo("`webshape` dtype")
                result = 0
            if self.webshape.ndim != 1:
                printInfo("`webshape` shape")
                result = 0
            if len(self.webshape) < 2:
                printInfo("`webshape` must has more than 2 layers")
                result = 0
        # inputdata
        if not isinstance(self.inputdata, ndarray):
            printInfo("`inputdata` type")
            result = 0
        else:
            if self.inputdata.dtype != float64:
                printInfo("`inputdata` dtype")
                result = 0
            if self.inputdata.ndim != 2:
                printInfo("`inputdata` shape")
                result = 0
        # outputdata
        if not isinstance(self.outputdata, ndarray):
            printInfo("`outputdata` type")
            result = 0
        else:
            if self.outputdata.dtype != float64:
                printInfo("`outputdata` dtype")
                result = 0
            if self.outputdata.ndim != 2:
                printInfo("`outputdata` shape")
                result = 0
        # trans
        if not isinstance(self.trans, ndarray):
            printInfo("`trans` type")
            result = 0
        else:
            if self.trans.dtype != object:
                printInfo("`trans` dtype")
                result = 0
            if self.trans.ndim != 1:
                printInfo("`trans` shape")
                result = 0
            for i in range(len(self.trans)):
                if not isinstance(self.trans[i], ndarray):
                    printInfo("`trans[{0}]` dtype".format(i))
                    result = 0
        # webshape
        if isinstance(self.inputdata, ndarray) and isinstance(self.webshape, ndarray):
            if self.inputdata.shape[1] != self.webshape[0]:
                printInfo("`inputdata` doesn't fit webshape[0]")
                result = 0
        if isinstance(self.outputdata, ndarray) and isinstance(self.webshape, ndarray):
            if self.outputdata.shape[1] != self.webshape[-1]:
                printInfo("`outputdata` doesn't fit webshape[-1]")
                result = 0
        if isinstance(self.inputdata, ndarray) and isinstance(self.outputdata, ndarray):
            if self.inputdata.shape[0] != self.outputdata.shape[0]:
                printInfo("`inputdata` doesn't fit `outputdata`")
                result = 0
        if isinstance(self.webshape, ndarray) and isinstance(self.trans, ndarray):
            for i in range(len(self.webshape) - 1):
                if self.trans[i].shape != (self.webshape[i] + 1, self.webshape[i + 1]):
                    printInfo("`trans`[{0}] shape wrong".format(i))
                    result = 0
        return bool(result)

# def keepsize(n, **kwargs):
#     '''return the sample size for a data set n.'''
#     return int(n - round((1 - 1. / n)**n * n))


# def split_train_test(a, keep=None, **kwargs):
#     '''Split the data to train set and test set by random sample.'''
#     n = a.shape[0]
#     if keep == None:
#         keep = keepsize(n)
#     keeplist = range(n)
#     removelist = []
#     while len(keeplist) > keep:
#         removelist.append(
#             keeplist.pop(np.random.random_integers(0, len(keeplist) - 1)))
#     np.random.shuffle(removelist)
#     np.random.shuffle(keeplist)
#     return a[keeplist], a[removelist]


# def toclass(a, classier=None, **kwargs):
#     '''convert a matrix to a matrix split every class.'''
#     if classier != None:
#         setitem = sorted(list(set([i.item() for i in a])))
#     else:
#         setitem = classier
#     result = np.matrix(np.zeros((a.shape[0], len(setitem))))
#     for i in range(len(setitem)):
#         result[:, i] = (a == setitem[i])
#     return result
