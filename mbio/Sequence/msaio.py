'''MSA and matrix IO.
'''

__author__ = 'Wenzhi Mao'
__all__ = ['ReadFasta', 'WriteFasta', 'ReadMatrix', 'WriteMatrix']

from numpy import dtype, array


def _Startup():
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


def ReadFasta(filename, **kwargs):
    '''This a function to read Fasta file.
    Given a filename and the function will return a list of sequences.'''
    from numpy import array
    from os.path import isfile
    from .msa import MSA
    from .Cfasta import readFasta
    if not isfile(filename):
        raise ValueError("File not found.")
    msa = readFasta(filename)
    msa = MSA(msa[1], labels=msa[0])
    return msa


def Fastaline(x, **kwargs):
    res = ''
    while len(x) != 0:
        res += ''.join(x[:60]) + '\n'
        x = x[60:]
    return res


def WriteFasta(filename, msa, **kwargs):
    '''This a function to write Fasta file.
    Given a filename and a MSA class to write.'''
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


def ReadMatrix(filename, dtype=None, l=None, **kwargs):
    '''Read Matrix if possible.'''
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
        if mlist[0].shape!=mlist[1].shape:
            raise ValueError("The shape is not the same.")
        return [mlist[0][0], mlist[1][0]]
    else:
        m = parsesinglematrix(filename)
        return m


def parsesinglematrix(filename, **kwargs):
    import struct
    from numpy import fromfile, isfile
    from mbio.Application import math as mmath
    from os.path import exists
    if not (exists(filename)):
    	raise ValueError("File {0} doesn't exists.".format(filename))
    if not isfile(filename):
    	raise ValueError("File {0} is not a file.".format(filename))
    tempmatrix = fromfile(filename)
    if int(tempmatrix.shape[0]**0.5)**2!=tempmatrix.shape[0]:
    	raise ValueError("Matrix from {0} is not a square.".format(filename))
    l=int(tempmatrix.shape[0]**0.5)
    tempmatrix.resize((l,l))
    return tempmatrix


def WriteMatrix(m, filename, **kwargs):
    try:
    	m.tofile(filename)
    except:
    	raise TypeError("Couldn't write Matrix.")

_Startup()
