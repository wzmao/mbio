'''MSA and matrix IO.
'''

__author__ = 'Wenzhi Mao'
__all__ = ['ReadFasta']

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

_Startup()
