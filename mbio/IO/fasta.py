'''Some Fasta file read and write functions.
'''

__author__ = 'Wenzhi Mao'

__all__ = ['ReadFasta']


def _Startup():
    '''Get _path__.'''
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


def ReadFasta(filename, **kwargs):
    '''This a function to read Fasta file.
    Given a filename and the function will return a list of sequences.'''
    from numpy import array
    from os.path import isfile
    from mbio.Sequence.msa import MSA
    if not isfile(filename):
        return ValueError("File not found.")
    from .c_fasta import readFasta
    msa = readFasta(filename)
    msa = MSA(msa[1], labels=msa[0])
    return msa


_Startup()
