# -*- coding: utf-8 -*-
"""This module contains the MRC file class.
"""

__author__ = 'Wenzhi Mao'

__all__ = ['MRC']


class MRCHeader():

    """A header class for mrc file."""

    def __init__(self, **kwargs):
        self.nx = self.ny = self.nz = None
        self.mode = None
        self.nxstart = self.nystart = self.nzstart = None
        self.mx = self.my = self.mz = None
        self.cella = [None] * 3
        self.cellb = [None] * 3
        self.mapc = None
        self.mapr = None
        self.maps = None
        self.dmin = self.dmax = self.dmean = None
        self.ispg = None
        self.nsymbt = None
        self.extra = None
        self.origin = [None] * 3
        self.map = None
        self.machst = None
        self.rms = None
        self.nlabels = None
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
                self.header = readHeader(filename=filename, mode=mode, header=self.header)
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
