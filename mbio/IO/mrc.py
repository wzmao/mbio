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
        self.xstart = self.ystart = self.zstart = None
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
                    from numpy import array, argsort
                    self = temp
                    for i in range(10):
                        self.label[i] = self.label[i][:80]
                        if self.label[i].find('\0') != -1:
                            self.label[i] = self.label[i][
                                :self.label[i].find("\0")]
                    if self.symdata:
                        self.symdata = self.symdata[:80]
                        if self.symdata.find('\0') != -1:
                            self.symdata = self.symdata[
                                :self.symdata.find('\0')]
                    if self.extra:
                        self.extra = self.extra[:80]
                        if self.extra.find('\0') != -1:
                            self.extra = self.extra[:self.extra.find('\0')]
                    self.xstart, self.ystart, self.zstart = array(
                        [self.nxstart, self.nystart, self.nzstart])[argsort([self.mapc, self.mapr, self.maps])]
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
                    for i in range(10):
                        self.label[i] = self.label[i][:80]
                        if self.label[i].find('\0') != -1:
                            self.label[i] = self.label[i][
                                :self.label[i].find("\0")]
                    if self.symdata:
                        self.symdata = self.symdata[:80]
                        if self.symdata.find('\0') != -1:
                            self.symdata = self.symdata[
                                :self.symdata.find('\0')]
                    if self.extra:
                        self.extra = self.extra[:80]
                        if self.extra.find('\0') != -1:
                            self.extra = self.extra[:self.extra.find('\0')]
            else:
                from .output import printError
                printError("The file doesn't exists or is not a file.")
        else:
            from .output import printError
            printError("The filename must be provided.")

    def printInfomation(self, **kwargs):
        """Print the information from the header."""

        from .output import printInfo as p

        p("Num of columns, rows and sections: {0} {1} {2}".format(
            self.nx, self.ny, self.nz))
        p("Mode:                              {0}".format(self.mode))
        p("Num of First column, row, section: {0} {1} {2}".format(
            self.nxstart, self.nystart, self.nzstart))
        p("Num of intervals along x, y, z:    {0} {1} {2}".format(
            self.mx, self.my, self.mz))
        p("Cell dimensions in angstroms:      {0:.2f} {1:.2f} {2:.2f}".format(
            self.cella[0], self.cella[1], self.cella[2]))
        p("Cell angles in degrees:            {0:.2f} {1:.2f} {2:.2f}".format(
            self.cellb[0], self.cellb[1], self.cellb[2]))
        p("Axis for cols, rows, sections:     {0} {1} {2}".format(
            self.mapc, self.mapr, self.maps))
        p("Min, max, mean density value:      {0:.6f} {1:.6f} {2:.6f}".format(
            self.dmin, self.dmax, self.dmean))
        p("Space group number:                {0}".format(self.ispg))
        p("Origin in X,Y,Z:                   {0:.4f} {1:.4f} {2:.4f}".format(
            self.origin[0], self.origin[1], self.origin[2]))
        p("Machine stamp:                     {0}".format(self.machst))
        p("rms deviationfrom mean density:    {0}".format(self.rms))
        p("Num of labels being used:          {0}".format(self.nlabels))
        if self.nlabels != 0:
            p("Labels:")
        for i in self.label:
            if i != "":
                p("\t{0}".format(i))
        p("Num of bytes for symmetry data:    {0}".format(self.nsymbt))
        if self.nsymbt != 0:
            p("\t{0}".format(self.symdata))

    def __repr__(self):

        return "MRCHeader"

    def setValue(self, label, value=None, **kwargs):
        """Set the value for a label."""

        setattr(self, label, value)

    def getValue(self, label, default=None, **kwargs):
        """Get the value for a label."""

        getattr(self, label, default)


class MRC():

    """This is a class to read and write MRC file.
    The data will always been store as x,y,z oreder."""

    def __init__(self, filename=None, **kwargs):
        """Parse data from the given file."""

        self.header = MRCHeader()
        self.data = None
        if filename:
            self.parseData(filename=filename)

    def __getattr__(self, name, **kwargs):

        if name in ['data', 'header']:
            return getattr(self, name)
        else:
            try:
                return getattr(self.header, name)
            except:
                return None

    def __setattr__(self, name, value, **kwargs):

        if name == 'data':
            self.__dict__[name] = value
        elif name == 'header':
            self.__dict__[name] = value
        else:
            if name in self.header.__dict__.keys():
                setattr(self.header, name, value)
            elif name in self.__dict__.keys():
                setattr(self, name, value)
            else:
                pass

    def __repr__(self):

        return "MRC"

    def __str__(self):

        return "MRC"

    def __dir__(self, **kwargs):

        return self.__dict__.keys() + self.header.__dict__.keys()

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
                from .output import printInfo, printError

                if getattr(self, 'header', None):
                    del self.header

                printInfo("Parsing the Header from file {0}.".format(filename))
                self.header = MRCHeader(filename=filename)

                if getattr(self, 'data', None):
                    printInfo("Some data exists already, overwrite it.")
                    del self.data

                if self.header.mode in [3, 4]:
                    printError(
                        "Sorry, we don't support the complex format yet.")
                    del self.data
                    self.data = None
                    return None
                else:
                    if self.header.mode == 0:
                        self.data = zeros(
                            (self.header.nz, self.header.ny, self.header.nx), dtype=int8)
                    elif self.header.mode == 1:
                        self.data = zeros(
                            (self.header.nz, self.header.ny, self.header.nx), dtype=int16)
                    elif self.header.mode == 2:
                        self.data = zeros(
                            (self.header.nz, self.header.ny, self.header.nx), dtype=float32)
                    elif self.header.mode == 5:
                        self.data = zeros(
                            (self.header.nz, self.header.ny, self.header.nx), dtype=uint8)
                    elif self.header.mode == 6:
                        self.data = zeros(
                            (self.header.nz, self.header.ny, self.header.nx), dtype=uint16)
                    else:
                        printError(
                            "Couldn't understand the mode {0}".format(self.header.mode))
                        del self.data
                        self.data = None
                        return None
                    printInfo(
                        "Parsing the Data from file {0}.".format(filename))
                    self.data = self.data - 1
                    temp = readData(
                        filename=filename, nsymbt=self.header.nsymbt,
                        datamode=self.header.mode, data=self.data,
                        size=self.header.nz * self.header.ny * self.header.nx)
                    if isinstance(temp, tuple):
                        del self.data
                        self.data = None
                        if temp[0] == None:
                            printError(temp[1])
                        else:
                            printError("Couldn't parse the Error information.")
                        return None
                    else:
                        from numpy import transpose, argsort
                        if set([self.header.mapc, self.header.mapr, self.header.maps]) != set([1, 2, 3]):
                            printError(
                                "The MRC header contains no clear axis.(mapc, mapr and maps must cotain all 1,2,3.)")
                            printError("Keep the data as it.")
                            self.data = temp
                            return None
                        else:
                            temporder = [
                                self.header.maps, self.header.mapr, self.header.mapc]
                            self.data = transpose(temp, argsort(temporder))
                            del temp
            else:
                printError("The file doesn't exists or is not a file.")
                return None
        else:
            printError("The filename must be provided.")
            return None

    def writeData(self, filename, skipupdate=False, force=False, **kwargs):
        """Write the MRC file into file.
        The header and data format will automaticly update.
        You could skip the update using `skipupdate` option.
        You could force it to overwrite files with `force` option."""

        from .output import printInfo, printError
        from os.path import exists, isfile
        from numpy import transpose, array

        if filename:
            if exists(filename):
                if not isfile(filename):
                    printError("The path is not a file.")
                    return None
                else:
                    if not force:
                        back = raw_input(
                            "* File {0} exists, do you want to overwrite it?(y/n)".format(filename))
                        while back.strip().lower() not in ['y', 'n']:
                            back = raw_input(
                                "* File {0} exists, do you want to overwrite it?(y/n)".format(filename))
                        if back.strip().lower() == 'n':
                            printInfo("File not write.")
                            return None
        else:
            printError("The filename must be provided.")
            return None
        if isinstance(self.data, type(None)):
            printError("No data to write.")
            return None
        find = False
        for i in range(10):
            if self.label[i].startswith("Written by mbio"):
                find = True
                from time import ctime
                from .. import __version__
                self.label[i] = "Written by mbio {0} {1}".format(
                    __version__, ctime())
                self.label = self.label[:i] + \
                    self.label[i + 1:] + [self.label[i]]
                self.label = [j for j in self.label if j != ""]
                self.label = self.label + [""] * (10 - len(self.label))
                break
        if not find:
            if self.nlabels != 10:
                from time import ctime
                from .. import __version__
                self.label[self.nlabels] = "Written by mbio {0} {1}".format(
                    __version__, ctime())
                self.nlabels += 1
        if not skipupdate:
            self.update()
        from .Cmrc import writeData
        if set([self.header.mapc, self.header.mapr, self.header.maps]) != set([1, 2, 3]):
            printError(
                "The MRC header contains no clear axis.(mapc, mapr and maps must cotain all 1,2,3.)")
            printError("Change it automaticly.")
            self.header.mapc, self.header.mapr, self.header.maps = 1, 2, 3
        self.header.nxstart, self.header.nystart, self.header.nzstart = array(
            [self.header.xstart, self.header.ystart, self.header.zstart])[[self.header.mapc - 1, self.header.mapr - 1, self.header.maps - 1]]
        printInfo("Writing MRC to {0}".format(filename))
        temp = writeData(header=self.header, data=transpose(
            self.data, (self.header.maps - 1, self.header.mapr - 1, self.header.mapc - 1)), filename=filename)
        if isinstance(temp, tuple):
            if temp[0] == None:
                printError(temp[1])
            else:
                printError("Couldn't parse the Error information.")
            return None
        elif temp == 0:
            return None
        else:
            printError("Couldn't parse the Error information.")

    def update(self, **kwargs):
        """Update the MRC header information from the data array.
        Update the MRC data format based on the `header.mode`
        Include: nx, ny, nz, dmin, dmax, dmean, rms, nsymbt, nlabels and sort label
            nxstart, nystart, nzstart.
        Correct mapc, mapr and maps automaticly."""

        from numpy import array, int8, int16, float32, uint8, uint16
        from .output import printError

        if set([self.header.mapc, self.header.mapr, self.header.maps]) != set([1, 2, 3]):
            printError(
                "The MRC header contains no clear axis.(mapc, mapr and maps must cotain all 1,2,3.)")
            printError("Change it automaticly.")
            self.header.mapc, self.header.mapr, self.header.maps = 1, 2, 3
        self.header.nx, self.header.ny, self.header.nz = array(
            self.data.shape)[[self.header.mapc - 1, self.header.mapr - 1, self.header.maps - 1]]
        self.header.nxstart, self.header.nystart, self.header.nzstart = array(
            [self.header.xstart, self.header.ystart, self.header.zstart])[[self.header.mapc - 1, self.header.mapr - 1, self.header.maps - 1]]
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

    def truncMatrix(self, index=[None,None,None,None,None,None] ,**kwargs):
        """Trunc the matrix by index. Related values will change accordingly.
        You need provide the start and end index(will be included) of x,y and z.
        Exapmle: 
            MRC.truncMatrix([xstart, xend, ystart, yend, zstart, zend])
        You could use *None* to indicate start from begin or to the end.
        """

        from .output import printError, printInfo
        from numpy import array
        if len(index)!=6:
            printError("Must provide 6 indeces.")
            return None
        if index==[None]*6:
            printInfo("Nothing changed.")
            return None
        xstart, xend, ystart, yend, zstart, zend = index
        if xstart==None:
            xstart=0
        if ystart==None:
            ystart=0
        if zstart==None:
            zstart=0
        if xend==None:
            xend=self.data.shape[0]+1
        else:
            xend+=1
        if yend==None:
            yend=self.data.shape[1]+1
        else:
            yend+=1
        if zend==None:
            zend=self.data.shape[2]+1
        else:
            zend+=1
        if not 0<=xstart<self.data.shape[0]:
            printError("xstart is not in the range of x.")
            return None
        if not 0<=xend<self.data.shape[0]:
            printError("xend is not in the range of x.")
            return None
        if not xstart<xend:
            printError("xstart must less than xend.")
            return None
        if not 0<=ystart<self.data.shape[0]:
            printError("ystart is not in the range of y.")
            return None
        if not 0<=yend<self.data.shape[0]:
            printError("yend is not in the range of y.")
            return None
        if not ystart<yend:
            printError("ystart must less than yend.")
            return None
        if not 0<=zstart<self.data.shape[0]:
            printError("zstart is not in the range of z.")
            return None
        if not 0<=zend<self.data.shape[0]:
            printError("zend is not in the range of z.")
            return None
        if not zstart<zend:
            printError("zstart must less than zend.")
            return None
        self.data=self.data[xstart:xend,ystart:yend,zstart:zend]
        self.header.xstart+=xstart
        self.header.ystart+=ystart
        self.header.zstart+=zstart
        self.header.nxstart, self.header.nystart, self.header.nzstart = array(
            [self.header.xstart, self.header.ystart, self.header.zstart])[[self.header.mapc - 1, self.header.mapr - 1, self.header.maps - 1]]

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

    def printInfomation(self, **kwargs):
        """Print the information from the header."""

        self.header.printInfomation()
