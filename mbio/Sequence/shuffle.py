# -*- coding: utf-8 -*-
"""This module contains some shuffle algorithms for protein sequence analysis.
Shuffle could be a good method to eliminate some false positive signals.
This module include the shuffle function for OMES, MI and MIp.
"""

__author__ = 'Wenzhi Mao'
__all__ = ['shuffleOMES', 'shuffleMI', 'shuffleMIp', 'shuffleAll']


def shuffleOMES(msa, times=10000, ambiguity=True, **kwargs):
    """Shuffle protein MSA within each column and calculate the p-value.
    Each OMES(Observed Minus Expected Squared) score will have a corresponding p-value.
    This function will return a numpy matrix with the p values.
    The core function is written in C."""

    from .correlation import getMSA
    msa = getMSA(msa)
    from .Cshuffle import shuffleomes
    length = msa.shape[1]
    from numpy import empty
    p = empty((length, length), dtype=float)
    p = shuffleomes(msa, p, ambiguity=bool(ambiguity), times=times)
    return p


def shuffleMI(msa, times=10000, ambiguity=True, **kwargs):
    """Shuffle protein MSA within each column and calculate the p-value.
    Each MI(Mutual Information) score will have a corresponding p-value.
    This function will return a numpy matrix with the p values.
    The core function is written in C."""

    from .correlation import getMSA
    msa = getMSA(msa)
    from .Cshuffle import shufflemi
    length = msa.shape[1]
    from numpy import empty
    p = empty((length, length), dtype=float)
    p = shufflemi(msa, p, ambiguity=bool(ambiguity), times=times)
    return p


def shuffleMIp(msa, times=10000, ambiguity=True, **kwargs):
    """Shuffle protein MSA within each column and calculate the p-value.
    Each MIp(Corrected Mutual Information) score will have a corresponding p-value.
    This function will return a numpy matrix with the p values.
    The core function is written in C."""

    from .correlation import getMSA
    msa = getMSA(msa)
    from .Cshuffle import shufflemip
    length = msa.shape[1]
    from numpy import empty
    p = empty((length, length), dtype=float)
    p = shufflemip(msa, p, ambiguity=bool(ambiguity), times=times)
    return p


def shuffleAll(msa, times=10000, ambiguity=True, mi=1, mip=1, omes=1, **kwargs):
    """Shuffle protein MSA within each column and calculate the p-value.
    This function combine MI, MIp and OMES.
    It will just shuffle once and return the p-value for all methods.
    This function will return a numpy matrix dict with the p values.
    The core function is written in C."""

    from .correlation import getMSA
    msa = getMSA(msa)
    from .Cshuffle import shuffleall
    length = msa.shape[1]
    mi = 1 if mi else 0
    mip = 1 if mip else 0
    omes = 1 if omes else 0
    from numpy import zeros, dtype
    pmi = zeros((length, length), dtype=float)
    pmip = zeros((length, length), dtype=float)
    pomes = zeros((length, length), dtype=float)
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
