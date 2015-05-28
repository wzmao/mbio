# -*- coding: utf-8 -*-
"""This module contains the MRC file class.
"""

__author__ = 'Wenzhi Mao'

__all__ = ['MRC']


class MRCHeader():

    """A header class for mrc file."""

    def __init__(self, filename=None, **kwargs):
        """Provide the filename to parse or set it later."""
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
        self.label = [None] * 10
        self.symdata = None
        if filename:
            from os.path import exists, isfile
            if exists(filename) and isfile(filename):
                from .Cmrc import readHeader
                temp = readHeader(filename=filename, header=self)
                if isinstance(temp, tuple):
                    from .output import printError
                    if temp[0] == None:
                        printError(temp[1])
                    else:
                        printError("Couldn't parse the Error information.")
                    return None
                else:
                    self = temp
            else:
                from .output import printError
                printError("The file doesn't exists or is not a file.")

    def parseHeader(self, filename=None, **kwargs):
        """Parse the MRC header information from the given file."""
        if filename:
            from os.path import exists, isfile
            if exists(filename) and isfile(filename):
                from .Cmrc import readHeader
                temp = readHeader(filename=filename, header=self)
                if isinstance(temp, tuple):
                    from .output import printError
                    if temp[0] == None:
                        printError(temp[1])
                    else:
                        printError("Couldn't parse the Error information.")
                    return None
                else:
                    self = temp
            else:
                from .output import printError
                printError("The file doesn't exists or is not a file.")
        else:
            from .output import printError
            printError("The filename must be provided.")

    def __repr__(self):
        return "MRCHeader"

    def setValue(self, label, value=None, **kwargs):
        """Set the value for a label."""
        setattr(self, label, value)

    def getValue(self, label, default=None, **kwargs):
        """Get the value for a label."""
        getattr(self, label, default)


class MRC():

    """This is a class to read and write MRC file."""

    def __init__(self, filename=None, **kwargs):
        """Parse data from the given file."""
        self.header = MRCHeader()
        self.data = None
        if filename:
            self.parseData(filename=filename)

    def parseHeader(self, filename=None, **kwargs):
        """Parse the header only from a given file.
        If the data will be parsed in the future, the header will be overwrited
        by the new data file's header."""
        if filename:
            from os.path import exists, isfile
            if exists(filename) and isfile(filename):
                self.header = MRCHeader(filename=filename)
            else:
                from .output import printError
                printError("The file doesn't exists or is not a file.")
        else:
            from .output import printError
            printError("The filename must be provided.")

    def parseData(self, filename=None, **kwargs):
        """Parse the data and header from a given file.
        If the header or data have already exists, all will be overwrited."""
        if filename:
            from os.path import exists, isfile
            if exists(filename) and isfile(filename):
                from .Cmrc import readData
                from numpy import zeros, int8, int16, float32, uint8, uint16
                from .output import printInfo

                if getattr(self, 'header', None):
                    del self.header

                printInfo("Parsing the Header from file {0}.".format(filename))
                self.header = MRCHeader(filename=filename)

                if getattr(self, 'data', None):
                    printInfo("Some data exists already, overwrite it.")
                    del self.data

                if self.header.mode in [3, 4]:
                    from .output import printError
                    printError(
                        "Sorry, we don't support the complex format yet.")
                    del self.data
                    self.data = None
                    return None
                else:
                    if self.header.mode == 0:
                        self.data = zeros(
                            (self.header.nx, self.header.ny, self.header.nz), dtype=int8)
                    elif self.header.mode == 1:
                        self.data = zeros(
                            (self.header.nx, self.header.ny, self.header.nz), dtype=int16)
                    elif self.header.mode == 2:
                        self.data = zeros(
                            (self.header.nx, self.header.ny, self.header.nz), dtype=float32)
                    elif self.header.mode == 5:
                        self.data = zeros(
                            (self.header.nx, self.header.ny, self.header.nz), dtype=uint8)
                    elif self.header.mode == 6:
                        self.data = zeros(
                            (self.header.nx, self.header.ny, self.header.nz), dtype=uint16)
                    else:
                        from .output import printError
                        printError(
                            "Couldn't understand the mode {0}".format(self.header.mode))
                        del self.data
                        self.data = None
                        return None
                    printInfo(
                        "Parsing the Data from file {0}.".format(filename))
                    self.data = self.data - 1.0
                    temp = readData(
                        filename=filename, nsymbt=self.header.nsymbt,
                        datamode=self.header.mode, data=self.data,
                        size=self.header.nx * self.header.ny * self.header.nz)
                    if isinstance(temp, tuple):
                        from .output import printError
                        del self.data
                        self.data = None
                        if temp[0] == None:
                            printError(temp[1])
                        else:
                            printError("Couldn't parse the Error information.")
                        return None
                    else:
                        self.data = temp
            else:
                from .output import printError
                printError("The file doesn't exists or is not a file.")
        else:
            from .output import printError
            printError("The filename must be provided.")

    def update(self, **kwargs):
        """Update the MRC header information from the data array.
        Update the MRC data format based on the `header.mode`
        Include: nx, ny, nz, dmin, dmax, dmean, rms, nsymbt, nlabels and sort label."""
        from numpy import array, int8, int16, float32, uint8, uint16
        self.header.nx, self.header.ny, self.header.nz = self.data.shape
        self.header.dmin = self.data.min()
        self.header.dmax = self.data.max()
        self.header.dmean = self.data.mean()
        self.header.rms = (((self.data - self.data.mean())**2).mean())**.5
        if self.header.symdata:
            self.header.nsymbt = 80
            self.header.symdata = self.header.symdata[:80]
        else:
            self.header.nsymbt = 0
            self.header.symdata = None
        self.header.nlabels = sum(
            [1 if i != "" else 0 for i in self.header.label])
        self.header.label = [i[:80] for i in self.header.label if i != ""]
        self.header.label = self.header.label + \
            [""] * (10 - len(self.header.label))
        if {0: int8, 1: int16, 2: float32, 5: uint8, 6: uint16}[self.header.mode] != self.data.dtype:
            self.data = array(self.data,
                              dtype={0: int8, 1: int16, 2: float32, 5: uint8, 6: uint16}[self.header.mode])

    def __getattr__(self, name, **kwargs):
        if name in ['data', 'header']:
            return getattr(self, name)
        else:
            try:
                return getattr(self.header, name)
            except:
                return None

    def __dir__(self, **kwargs):
        return self.__dict__.keys() + self.header.__dict__.keys()

    def getMatrixShape(self, **kwargs):
        """Get the data shape from the header information.
        Caution: it could be different with the data array."""
        if (isinstance(self.header.nx, int) and
                isinstance(self.header.ny, int) and isinstance(self.header.nz, int)):
            return (self.header.nx, self.header.ny, self.header.nz)
        else:
            from .output import printError
            printError("There is no header information here.")
            return None

    def getArray(self, **kwargs):
        """Get the data from the MRC class"""
        return self.data

    def setMode(self, mode=2, **kwargs):
        """Set the data format for the data.
        The data will be change the format accordingly.
        Data type :
            0       image : signed 8-bit bytes range -128 to 127
            1       image : 16-bit halfwords
            2       image : 32-bit reals
            3       transform : complex 16-bit integers (not support now)
            4       transform : complex 32-bit reals (not support now)
            5       image : unsigned 8-bit range 0 to 255
            6       image : unsigned 16-bit range 0 to 65535"""
        from numpy import array, int8, int16, float32, uint8, uint16
        from .output import printError
        if mode not in range(7):
            printError("Mode must be 0,1,2,3,4,5,6.")
        elif mode in [3, 4]:
            printError("Sorry, the complex format is not supported now.")
        self.header.mode = mode
        if {0: int8, 1: int16, 2: float32, 5: uint8, 6: uint16}[self.header.mode] != self.data.dtype:
            self.data = array(self.data,
                              dtype={0: int8, 1: int16, 2: float32, 5: uint8, 6: uint16}[self.header.mode])
