'''Shuffle could be a good method to eliminate some false positive signals.
This module include the shuffle function for OMES, MI and MIp.
'''

__author__ = 'Wenzhi Mao'
__all__ = ['ShuffleOMES', 'ShuffleMI', 'ShuffleMIp', 'ShuffleAll']

from .correlation import getMSA
from numpy import dtype, zeros, empty, ones


def _Startup():
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


def ShuffleOMES(msa, times=10000, ambiguity=True, **kwargs):
    msa = getMSA(msa)
    from .Cshuffle import shuffleomes
    length = msa.shape[1]
    p = empty((length, length), float)
    p = shuffleomes(msa, p, ambiguity=bool(ambiguity), times=times)
    return p


def ShuffleMI(msa, times=10000, ambiguity=True, **kwargs):
    msa = getMSA(msa)
    from .Cshuffle import shufflemi
    length = msa.shape[1]
    p = empty((length, length), float)
    p = shufflemi(msa, p, ambiguity=bool(ambiguity), times=times)
    return p


def ShuffleMIp(msa, times=10000, ambiguity=True, **kwargs):
    msa = getMSA(msa)
    from .Cshuffle import shufflemip
    length = msa.shape[1]
    p = empty((length, length), float)
    p = shufflemip(msa, p, ambiguity=bool(ambiguity), times=times)
    return p


def ShuffleAll(msa, times=10000, ambiguity=True, mi=1, mip=1, omes=1, **kwargs):
    msa = getMSA(msa)
    from .Cshuffle import shuffleall
    length = msa.shape[1]
    mi = 1 if mi else 0
    mip = 1 if mip else 0
    omes = 1 if omes else 0
    pmi = zeros((length, length), float)
    pmip = zeros((length, length), float)
    pomes = zeros((length, length), float)
    p = shuffleall(msa, ambiguity=bool(ambiguity), times=times,
                   hasmi=mi, hasmip=mip, hasomes=omes,
                   pmi=pmi, pmip=pmip, pomes=pomes)
    if not mi:
        p.pop('MI')
    if not mip:
        p.pop('MIp')
    if not omes:
        p.pop('OMES')
    return p

_Startup()
