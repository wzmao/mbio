# -*- coding: utf-8 -*-
"""This module contains some protein sequence correlation algorithms.
Protein sequence correlation analysis could provide the evolutionary 
information about structure, contact and so on.
Now, this module includes MI, MIp, OMES, SCA and DI.
There are also APC and BND correction to refine the correlation matrix further.
"""

__author__ = 'Wenzhi Mao'
__all__ = ['buildMI', 'buildMIp', 'buildOMES',
           'buildSCA', 'buildDI', 'buildDCA',
           'calcMeff', 'applyAPC', 'applyBND']


def getMSA(msa):
    """Return MSA character array."""

    from numpy import dtype
    try:
        msa = msa._getArray()
    except AttributeError:
        pass
    try:
        msa = msa.seq
    except AttributeError:
        pass
    try:
        dtype_, ndim, shape = msa.dtype, msa.ndim, msa.shape
    except AttributeError:
        raise TypeError('msa must be an MSA instance or a 2D character array')
    if dtype_ != dtype('|S1') or ndim != 2:
        raise TypeError('msa must be an MSA instance or a 2D character array')
    return msa


def buildMI(msa, ambiguity=True, turbo=True, **kwargs):
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
    :func:`applyAPC` and :func:`applyBND` methods, respectively.

    [GGB05] Gloor, Gregory B., et al. "Mutual information in protein multiple 
    sequence alignments reveals two classes of coevolving positions." 
    Biochemistry 44.19 (2005): 7156-7165."""

    msa = getMSA(msa)

    from numpy import empty
    from .Ccorrelation import msamutinfo
    length = msa.shape[1]
    mutinfo = empty((length, length), float)
    mutinfo = msamutinfo(msa, mutinfo,
                         ambiguity=bool(ambiguity), turbo=bool(turbo),
                         norm=bool(kwargs.get('norm', False)),
                         debug=bool(kwargs.get('debug', False)))

    return mutinfo


def buildMIp(msa, ambiguity=True, turbo=True, **kwargs):
    """It is a function to calculate the MIp matrix.
    MIp is a APC corrected form of MI matrix.
    It could provide better information than MI.

    [DSD08] Dunn SD, Wahl LM, Gloor GB. Mutual information without the 
    influence of phylogeny or entropy dramatically improves residue 
    contact prediction. Bioinformatics 2008 24(3):333-340."""

    return ApplyAPC(buildMI(msa, ambiguity=True, turbo=True, **kwargs))


def buildOMES(msa, ambiguity=True, turbo=True, **kwargs):
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
    are considered as gaps.

    [KI02] Kass, Itamar, and Amnon Horovitz. "Mapping pathways of allosteric 
    communication in GroEL by analysis of correlated mutations." Proteins: 
    Structure, Function, and Bioinformatics 48.4 (2002): 611-617."""

    msa = getMSA(msa)

    from numpy import empty
    from .Ccorrelation import msaomes
    length = msa.shape[1]
    omes = empty((length, length), float)
    omes = msaomes(msa, omes, ambiguity=bool(ambiguity), turbo=bool(turbo),
                   debug=bool(kwargs.get('debug', False)))

    return omes


def buildSCA(msa, turbo=True, **kwargs):
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
    are considered as gaps.

    [LSW99] Lockless, Steve W., and Rama Ranganathan. "Evolutionarily conserved 
    pathways of energetic connectivity in protein families." Science 286.5438 
    (1999): 295-299.

    [HN09] Halabi, Najeeb, et al. "Protein sectors: evolutionary units of 
    three-dimensional structure." Cell 138.4 (2009): 774-786."""

    msa = getMSA(msa)

    from numpy import zeros
    from .Ccorrelation import msasca
    length = msa.shape[1]
    sca = zeros((length, length), float)
    sca = msasca(msa, sca, turbo=bool(turbo))
    return sca


def buildDI(msa, seqid=.8, pseudo_weight=.5, refine=False,
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

    [WM09] Weigt, Martin, et al. "Identification of direct residue contacts 
    in protein–protein interaction by message passing." Proceedings of the 
    National Academy of Sciences 106.1 (2009): 67-72.

    [MF11] Morcos, Faruck, et al. "Direct-coupling analysis of residue 
    coevolution captures native contacts across many protein families." 
    Proceedings of the National Academy of Sciences 108.49 
    (2011): E1293-E1301."""

    msa = getMSA(msa)

    from numpy import zeros
    from .Ccorrelation import msadipretest, msadirectinfo1, msadirectinfo2
    from numpy import matrix

    refine = 1 if refine else 0
    # msadipretest get some parameter from msa to set matrix size
    length, q = msadipretest(msa, refine=refine)
    c = matrix.dot(matrix(zeros((length * q, 1), float)),
                   matrix(zeros((1, length * q), float)))
    prob = zeros((length, q + 1), float)
    # msadirectinfo1 return c to be inversed and prob to be used
    meff, n, length, c, prob = msadirectinfo1(msa, c, prob, theta=1. - seqid,
                                              pseudocount_weight=pseudo_weight,
                                              refine=refine, q=q + 1)

    c = c.I

    di = zeros((length, length), float)
    # get final DI
    di = msadirectinfo2(n, length, c, prob, di, q + 1)
    del prob, c
    return di


def buildDCA(msa, seqid=.8, pseudo_weight=.5, refine=False,
             **kwargs):
    """Same as buildDI."""

    return buildDI(msa, seqid, pseudo_weight, refine, **kwargs)


def calcMeff(msa, seqid=.8, refine=False, weight=False, **kwargs):
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

    The weight for each sequence are returned when *weight* is **True**.

    [WM09] Weigt, Martin, et al. "Identification of direct residue contacts 
    in protein–protein interaction by message passing." Proceedings of the 
    National Academy of Sciences 106.1 (2009): 67-72.

    [MF11] Morcos, Faruck, et al. "Direct-coupling analysis of residue 
    coevolution captures native contacts across many protein families." 
    Proceedings of the National Academy of Sciences 108.49 
    (2011): E1293-E1301."""

    msa = getMSA(msa)

    from numpy import zeros
    from .Ccorrelation import msameff
    refine = 1 if refine else 0
    weight = 0 if weight else 1  # A Mark for return weighted array.
    if (not weight):
        w = zeros((msa.shape[0]), float)
        meff = msameff(msa, theta=1. - seqid, meff_only=weight,
                       refine=refine, w=w)
    else:
        meff = msameff(msa, theta=1. - seqid, meff_only=weight, refine=refine)
    return meff


def applyAPC(mutinfo, **kwargs):
    """Return a copy of *mutinfo* array after average product correction
    (default) or average sum correction is applied.

    [DSD08] Dunn SD, Wahl LM, Gloor GB. Mutual information without the 
    influence of phylogeny or entropy dramatically improves residue 
    contact prediction. Bioinformatics 2008 24(3):333-340."""

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


def applyBND(mat, **kwargs):
    """Return a BND refinement of a symetric matrix.
    If the matrix is not symetric. The calculation for mat+mat.T will be 
    performanced.
    It has comparable speed with the original MATLAB code.

    [SUP15] Sun, Hai‐Ping, et al. "Improving accuracy of protein contact 
    prediction using balanced network deconvolution." Proteins: Structure, 
    Function, and Bioinformatics (2015)."""

    try:
        ndim, shape = mat.ndim, mat.shape
    except AttributeError:
        raise TypeError('Matrix must be a 2D square array')

    from numpy import fill_diagonal
    from scipy.linalg import eigh
    from scipy.linalg.blas import dgemm
    if ndim != 2 or shape[0] != shape[1]:
        raise ValueError('Matrix must be a 2D square array')

    if mat.min() != mat.max():
        mat = (mat - mat.min()) / (mat.max() - mat.min())
    else:
        from ..IO.output import printInfo
        printInfo("The input matrix is a constant matrix.")
        return mat
    n = mat.shape[0]
    fill_diagonal(mat, 0.)
    mat_th = (mat + mat.T) / 2.
    d, u = eigh(mat_th)
    for i in range(n):
        if d[i] != 0:
            d[i] = (-1. + (1 + 4 * d[i] * d[i])**.5) / 2. / d[i]
    mat_new1 = dgemm(alpha=1.0, a=(u * d), b=u, trans_b=True)
    # mat_new1 = (u*d).dot(u.T) #old numpy dot,slower
    ind_edges = (mat_th > 0) * 1.0
    ind_nonedges = (mat_th == 0) * 1.0
    m1 = (mat * ind_nonedges).max()
    m2 = mat_new1.min()
    mat_new2 = (mat_new1 + max(m1 - m2, 0)) * ind_edges + (mat * ind_nonedges)
    m1 = mat_new2.min()
    m2 = mat_new2.max()
    if m1 != m2:
        mat_bnd = (mat_new2 - m1) / (m2 - m1)
    else:
        mat_bnd = mat_new2
    return mat_bnd
