# -*- coding: utf-8 -*-
"""This module contains the MRC file class.
"""

__author__ = 'Wenzhi Mao'

__all__ = ['MRC']


class MRCHeader():

    """A header class for mrc file."""

    def __init__(self, **kwargs):
        self.nx = self.ny = self.nz = 0
        self.mode = 2
        self.nxstart = self.nystart = self.nzstart = 0
        self.mx = self.my = self.mz = 0
        self.cella = [0.] * 3
        self.cellb = [90.] * 3
        self.mapc = 1
        self.mapr = 2
        self.maps = 3
        self.dmin = self.dmax = self.dmean = 0.
        self.ispg = 0
        self.nsymbt = 0
        self.extra = "\0"
        self.origin = [0.] * 3
        self.maps = "MAP\0"
        self.machst = 0
        self.rms = 0.
        self.nlabels = 0
        self.labels = [""] * 10

    def __repr__(self):
        return "MRCHeader"

    def setValue(self, label, value=None, **kwargs):
        setattr(self, label, value)

    def getValue(self, label, default=None, **kwargs):
        getattr(self, label, default)


class MRC():

    """This is a class to read and write MRC file."""

    def __init__(self, filename=None, **kwargs):
        if filename:
            self.parseHeader(filename)
            # self.parseData(filename)

    def parseHeader(self, filename=None, mode="r", **kwargs):
        if filename:
            from os.path import exists, isfile
            if exists(filename) and isfile(filename):
                from .Cmrc import readHeader
                self.header = MRCHeader()
                readHeader(filename=filename, mode=mode, header=self.header)
            else:
                from .output import printError
                printError("The file doesn't exists or is not a file.")
        else:
            from .output import printError
            printError("The filename must be provided.")

    # def parseData(self,filename=None,**kwargs):
    # 	if filename:
    # 		from os.path import exists, isfile
    # 		if exists(filename) and isfile(filename):
    # 			from .Cmrc import readData
    # 			self.header=MRCHeader()
    # 		else:
           #  		from .output import printError
            # printError("The file doesn't exists or is not a file.")
    # 	else:
    # 		from .output import printError
    # 		printError("The filename must be provided.")
