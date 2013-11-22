'''Different algorithm could identify protein position correlation.
We have MI, MIp, OMES and SCA now.
'''

__author__ = 'Wenzhi Mao'
__all__ = ['CalcMI', 'CalcMIp', 'CalcOMES', 'CalcSCA', 'CalcDI', 'CalcMeff', 'apc']

from numpy import dtype, zeros, empty, ones

def _Startup():
    '''Get _path__.'''
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


def getMSA(msa):
    """Return MSA character array."""

    try:
        msa = msa._getArray()
    except AttributeError:
        pass
    try:
        dtype_, ndim, shape = msa.dtype, msa.ndim, msa.shape
    except AttributeError:
        raise TypeError('msa must be an MSA instance or a 2D character array')
    if dtype_ != dtype('|S1') or ndim != 2:
        raise TypeError('msa must be an MSA instance or a 2D character array')
    return msa


def CalcMI(msa, ambiguity=True, turbo=True, **kwargs):
    """Return mutual information matrix calculated for *msa*, which may be an
    :class:`.MSA` instance or a 2D Numpy character array.  Implementation
    is case insensitive and handles ambiguous amino acids as follows:

      * **B** (Asx) count is allocated to *D* (Asp) and *N* (Asn)
      * **Z** (Glx) count is allocated to *E* (Glu) and *Q* (Gln)
      * **J** (Xle) count is allocated to *I* (Ile) and *L* (Leu)
      * **X** (Xaa) count is allocated to the twenty standard amino acids
      * Joint probability of observing a pair of ambiguous amino acids is
        allocated to all potential combinations, e.g. probability of **XX**
        is allocated to 400 combinations of standard amino acids, similarly
        probability of **XB** is allocated to 40 combinations of *D* and *N*
        with the standard amino acids.

    Selenocysteine (**U**, Sec) and pyrrolysine (**O**, Pyl) are considered
    as distinct amino acids.  When *ambiguity* is set **False**, all alphabet
    characters as considered as distinct types.  All non-alphabet characters
    are considered as gaps.

    Mutual information matrix can be normalized or corrected using
    :func:`applyMINormalization` and :func:`applyMICorrection` methods,
    respectively.  Normalization by joint entropy can performed using this
    function with *norm* option set **True**."""

    msa = getMSA(msa)

    from .c_correlation import msamutinfo
    length = msa.shape[1]
    mutinfo = empty((length, length), float)
    mutinfo = msamutinfo(msa, mutinfo,
                         ambiguity=bool(ambiguity), turbo=bool(turbo),
                         norm=bool(kwargs.get('norm', False)),
                         debug=bool(kwargs.get('debug', False)))

    return mutinfo


def CalcMIp(msa, ambiguity=True, turbo=True, **kwargs):
    '''It is a function to calculate the MIp matrix based on language C.
    Given the sequences in a list with no format.
    '''
    return apc(CalcMI(msa, ambiguity=True, turbo=True, **kwargs))


def CalcOMES(msa, ambiguity=True, turbo=True, **kwargs):
    """Return OMES (Observed Minus Expected Squared) covariance matrix
    calculated for *msa*, which may be an :class:`.MSA` instance or a 2D
    NumPy character array. OMES is defined as::

                        (N_OBS - N_EX)^2              (f_i,j - f_i * f_j)^2
      OMES_(i,j) = sum(------------------) = N * sum(-----------------------)
                             N_EX                           f_i * f_j

    Implementation is case insensitive and handles ambiguous amino acids
    as follows:

      * **B** (Asx) count is allocated to *D* (Asp) and *N* (Asn)
      * **Z** (Glx) count is allocated to *E* (Glu) and *Q* (Gln)
      * **J** (Xle) count is allocated to *I* (Ile) and *L* (Leu)
      * **X** (Xaa) count is allocated to the twenty standard amino acids
      * Joint probability of observing a pair of ambiguous amino acids is
        allocated to all potential combinations, e.g. probability of **XX**
        is allocated to 400 combinations of standard amino acids, similarly
        probability of **XB** is allocated to 40 combinations of *D* and *N*
        with the standard amino acids.

    Selenocysteine (**U**, Sec) and pyrrolysine (**O**, Pyl) are considered
    as distinct amino acids.  When *ambiguity* is set **False**, all alphabet
    characters as considered as distinct types.  All non-alphabet characters
    are considered as gaps."""

    msa = getMSA(msa)

    from .c_correlation import msaomes
    length = msa.shape[1]
    omes = empty((length, length), float)
    omes = msaomes(msa, omes, ambiguity=bool(ambiguity), turbo=bool(turbo),
                   debug=bool(kwargs.get('debug', False)))

    return omes

def CalcSCA(msa, turbo=True, **kwargs):
    """Return SCA matrix calculated for *msa*, which may be an :class:`.MSA`
    instance or a 2D Numpy character array.

    Implementation is case insensitive and handles ambiguous amino acids
    as follows:

      * **B** (Asx) count is allocated to *D* (Asp) and *N* (Asn)
      * **Z** (Glx) count is allocated to *E* (Glu) and *Q* (Gln)
      * **J** (Xle) count is allocated to *I* (Ile) and *L* (Leu)
      * **X** (Xaa) count is allocated to the twenty standard amino acids
      * Joint probability of observing a pair of ambiguous amino acids is
        allocated to all potential combinations, e.g. probability of **XX**
        is allocated to 400 combinations of standard amino acids, similarly
        probability of **XB** is allocated to 40 combinations of *D* and *N*
        with the standard amino acids.

    Selenocysteine (**U**, Sec) and pyrrolysine (**O**, Pyl) are considered
    as distinct amino acids.  When *ambiguity* is set **False**, all alphabet
    characters as considered as distinct types.  All non-alphabet characters
    are considered as gaps."""

    msa = getMSA(msa)
    from .c_correlation import msasca
    length = msa.shape[1]
    sca = zeros((length, length), float)
    sca = msasca(msa, sca, turbo=bool(turbo))
    return sca

def CalcDI(msa, seqid=.8, pseudo_weight=.5, refine=False,
                          **kwargs):
    """Return direct information matrix calculated for *msa*, which may be an
    :class:`.MSA` instance or a 2D Numpy character array.

    Sequences sharing sequence identity of *seqid* or more with another
    sequence are regarded as similar sequences for calculating their weights
    using :func:`.calcMeff`.

    *pseudo_weight* are the weight for pseudo count probability.

    Sequences are not refined by default. When *refine* is set **True**,
    the MSA will be refined by the first sequence and the shape of direct
    information matrix will be smaller.
    """

    msa = getMSA(msa)
    from .c_correlation import msadipretest, msadirectinfo1, msadirectinfo2
    from numpy import matrix

    refine = 1 if refine else 0
    # msadipretest get some parameter from msa to set matrix size
    length, q = msadipretest(msa, refine=refine)
    c = matrix.dot(matrix(zeros((length*q, 1), float)),
                   matrix(zeros((1, length*q), float)))
    prob = zeros((length, q+1), float)
    # msadirectinfo1 return c to be inversed and prob to be used
    meff, n, length, c, prob = msadirectinfo1(msa, c, prob, theta=1.-seqid,
                                              pseudocount_weight=pseudo_weight,
                                              refine=refine, q=q+1)

    c = c.I

    di = zeros((length, length), float)
    # get final DI
    di = msadirectinfo2(n, length, c, prob, di, q+1)
    del prob, c
    return di


def CalcMeff(msa, seqid=.8, refine=False, weight=False, **kwargs):
    """Return the Meff for *msa*, which may be an :class:`.MSA`
    instance or a 2D Numpy character array.

    Since similar sequences in an *msa* decreases the diversity of *msa*,
    *Meff* gives a weight for sequences in the *msa*.

    For example: One sequence in MSA has 5 other similar sequences in this
    MSA(itself included). The weight of this sequence is defined as 1/5=0.2.
    Meff is the sum of all sequence weights. In another word, Meff can be
    understood as the effective number of independent sequences.

    Sequences sharing sequence identity of *seqid* or more with another
    sequence are regarded as similar sequences to calculate Meff.

    Sequences are not refined by default. When *refine* is set **True**, the
    MSA will be refined by the first sequence.

    The weight for each sequence are returned when *weight* is **True**."""

    msa = getMSA(msa)
    from .c_correlation import msameff
    refine = 1 if refine else 0
    weight = 0 if weight else 1  # A Mark for return weighted array.
    if (not weight):
        w = zeros((msa.shape[0]), float)
        meff = msameff(msa, theta=1.-seqid, meff_only=weight,
                       refine=refine, w=w)
    else:
        meff = msameff(msa, theta=1.-seqid, meff_only=weight, refine=refine)
    return meff


def apc(mutinfo):
    """Return a copy of *mutinfo* array after average product correction
    (default) or average sum correction is applied."""

    try:
        ndim, shape = mutinfo.ndim, mutinfo.shape
    except AttributeError:
        raise TypeError('mutinfo must be a 2D square array')

    if ndim != 2 or shape[0] != shape[1]:
        raise ValueError('mutinfo must be a 2D square array')

    avg_mipos = mutinfo.sum(1) / (shape[0] - 1)
    avg_mi = avg_mipos.mean()

    mi = mutinfo.copy()
    for i, i_avg in enumerate(avg_mipos):
        for j, j_avg in enumerate(avg_mipos):
            mi[i, j] -= (i_avg * j_avg) / avg_mi
    return mi

_Startup()
