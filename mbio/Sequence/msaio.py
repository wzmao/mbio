"""This module contains some MSA and matrix IO functions.
"""

__author__ = 'Wenzhi Mao'
__all__ = ['readFasta', 'writeFasta', 'readMatrix', 'writeMatrix']


def readFasta(filename, **kwargs):
    """Read FASTA file.
    Given a filename and the function will return a :class:`MSA` of sequences.
    *** In plan: read from zipped files. ***"""

    from numpy import array
    from os.path import isfile
    from .msa import MSA
    from .Cfasta import readFasta
    if not isfile(filename):
        from ..IO.output import printError
        printError("File not found.")
        raise ValueError("File not found.")
    msa = readFasta(filename)
    msa = MSA(msa[1], labels=msa[0])
    return msa


def Fastaline(x, **kwargs):
    """Convert a sequence to FASTA format(60 per line)"""

    res = ''
    while len(x) != 0:
        res += ''.join(x[:60]) + '\n'
        x = x[60:]
    return res


def writeFasta(filename, msa, **kwargs):
    """Write FASTA file.
    Given a filename and a :class:`MSA` to write.
    *** In plan: write to zipped files. ***"""

    from os.path import exists
    if exists(filename):
        p = 1
        while (exists(filename + '({0})'.format(p))):
            p += 1
        filename = filename + '({0})'.format(p)
    f = open(filename, 'w')
    for i in xrange(msa.numseq):
        f.write('>{0}\n'.format(msa.label[i]))
        f.write('{0}'.format(Fastaline(msa.seq[i])))
    f.close()
    from ..IO.output import printInfo
    printInfo("The file {0} has been saved.".format(filename))


def readMatrix(filename, dtype=None, l=None, **kwargs):
    """Read matrix or matrices from file.
    Return a numpy array or a list of numpy array."""

    if (isinstance(filename, (tuple, list))):
        if len(filename) == 1:
            filename = filename[0]
        elif len(filename) > 1:
            filename = list(filename)[:2]
        else:
            filename = None
    elif type(filename) == str:
        filename = filename
    else:
        filename = None
    if not filename:
        raise ValueError("Couldn't understand the filename.")
    elif type(filename) == list:
        mlist = [parsesinglematrix(i) for i in filename]
        if mlist[0].shape != mlist[1].shape:
            raise ValueError("The shape is not the same.")
        return mlist
    else:
        m = parsesinglematrix(filename)
        return m


def parsesinglematrix(filename, **kwargs):
    """Parse matrix from file. Check if the matrix is square."""

    from numpy import fromfile
    from ..Application.math import issquare
    from os.path import exists, isfile
    if not (exists(filename)):
        raise ValueError("File {0} doesn't exists.".format(filename))
    if not isfile(filename):
        raise ValueError("File {0} is not a file.".format(filename))
    tempmatrix = fromfile(filename)
    if not issquare(tempmatrix.shape[0]):
        from ..IO.output import printError
        printError("Matrix from {0} is not a square.".format(filename))
        raise ValueError("Matrix from {0} is not a square.".format(filename))
    l = int(tempmatrix.shape[0] ** 0.5)
    tempmatrix.resize((l, l))
    return tempmatrix


def writeMatrix(m, filename, **kwargs):
    """Write a matrix to binary format using the numpy :func:`tofile`"""

    try:
        m.tofile(filename)
        from ..IO.output import printInfo
        printInfo("The file {0} has been saved.".format(filename))
    except:
        from ..IO.output import printError
        printError("Couldn't write Matrix.")
        raise TypeError("Couldn't write Matrix.")
