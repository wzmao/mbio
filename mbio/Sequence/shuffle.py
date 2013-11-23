'''Shuffle could be a good method to eliminate some false positive signals.
This module include the shuffle function for OMES, MI and MIp.
'''

__author__ = 'Wenzhi Mao'
__all__ = ['ShuffleOMES', 'ShuffleMI', 'ShuffleMIp']

from mbio.Sequence.correlation import getMSA
from numpy import dtype, zeros, empty, ones

def _Startup():
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


def ShuffleOMES(msa, times=10000, ambiguity=True):
    msa = getMSA(msa)
    from .c_shuffle import shuffleomes
    length = msa.shape[1]
    p = empty((length, length), float)
    p = shuffleomes(msa, p, ambiguity=bool(ambiguity), times=times)
    return p


def ShuffleMI(msa, times=10000, ambiguity=True):
    msa = getMSA(msa)
    from .c_shuffle import shufflemi
    length = msa.shape[1]
    p = empty((length, length), float)
    p = shufflemi(msa, p, ambiguity=bool(ambiguity), times=times)
    return p


def ShuffleMIp(msa, times=10000, ambiguity=True):
    msa = getMSA(msa)
    from .c_shuffle import shufflemip
    length = msa.shape[1]
    p = empty((length, length), float)
    p = shufflemip(msa, p, ambiguity=bool(ambiguity), times=times)
    return p


_Startup()
