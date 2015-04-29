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
           'calcMeff', 'applyAPC', 'applyBND',
           'buildPSICOV', 'applyPPV']


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
    # from scipy.linalg.lapack import dgetrf, dgetri
    # from scipy.linalg.lapack import dpotri

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

    # Another way:faster than normal,slower than atlas
    # d, e = dgetrf(c)[:2]
    # c = dgetri(d, e)[0]

    # c=spotri(c,)[0]

    di = zeros((length, length), float)
    # get final DI
    di = msadirectinfo2(n, length, c, prob, di, q + 1)
    del prob, c

    return di


def buildDCA(msa, seqid=.8, pseudo_weight=.5, refine=False,
             **kwargs):
    """The DCA matrix function see `buildDI` for more information."""

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
    from scipy.linalg.lapack import dsyevr
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
    # Double symetric eigvector relatively robust representation (RRR)
    [d, u, t] = dsyevr(mat_th)
    for i in range(n):
        if d[i] != 0:
            d[i] = (-1. + (1 + 4 * d[i] * d[i])**.5) / 2. / d[i]
    # Double general matrix multiply
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


def buildPSICOV_expert(msa,
                       approx_lasso_flag=0,
                       pre_shrink_flag=1,
                       force=0,
                       use_raw_not_ppv=1,
                       apply_apc_flg=1,
                       set_default_rho=-1,
                       target_fraction_of_none_zero=0,
                       set_lasso_convergence_threshold=1e-4,
                       identity_threshold=-1,
                       set_pseudocount_value=1,
                       set_minimum_sequence_separation=5,
                       rho_parameter_file="",
                       maximum_fraction_of_gaps=0.9,
                       maxthread=-1,
                       **kwargs):
    """Return PSICOV matrix calculated for *msa*, which may be an
    :class:`.MSA` instance or a 2D Numpy character array.

    You could use `buildPSICOV` for normal usage. `buildPSICOV_expert` provide
    more options.

    Options:
        `approx_lasso_flag`: use approximate LASSO algorithm, default as 0.
        `pre_shrink_flag`: Shrink the sample covariance matrix towards 
            shrinkage target F = Diag(1,1,1,...,1), default as 1.
        `force`: force the calculation if the sequence is not sufficient,
            default as 0.
        `use_raw_not_ppv`: return the raw PSICOV result instead of the PPV
            result(A model from [JDT12]), default as 1.
        `apply_apc_flg`: apply the APC after the calculation, default as 1.
            The result after APC could be negtive, but the original lowest
            score remains 0.
        `set_default_rho`: set the initial rho parameter, default as -1 (not
            specified)
        `target_fraction_of_none_zero`: set the target matrix density value 
            (none-zero fraction), should be in range 5e-5 - 1. default as 0
            (not specified)
        `set_lasso_convergence_threshold`: set Lasso convergence threshold,
            default as 1e-4.
        `identity_threshold`: select BLOSUM-like weighting with given identity
            threshold, default as -1(selects threshold automatically)
        `set_pseudocount_value`: set pseudocount value, default as 1.
        `set_minimum_sequence_separation`: set pseudocount value, default as 5.
        `rho_parameter_file`: give the rho parameter file, default ""(no file)
        `maximum_fraction_of_gaps`: maximum fraction of gaps, default as 0.9.
        `maxthread`: the max number of threads, default as -1(choose max).

    [JDT12] Jones, David T., et al. "PSICOV: precise structural contact
    prediction using sparse inverse covariance estimation on large multiple
    sequence alignments." Bioinformatics 28.2 (2012): 184-190."""

    msa = getMSA(msa)

    if (not 5e-5 < float(target_fraction_of_none_zero) < 1) and (target_fraction_of_none_zero != 0):
        from ..IO.output import printError
        printError("target_fraction_of_none_zero must between 5e-5 and 1.")
        return None
    from numpy import zeros
    from .Ccorrelation_p import msapsicov
    length = msa.shape[1]
    psicov = zeros((length, length), float)
    psicov = msapsicov(msa, psicov,
                       approxflg=bool(approx_lasso_flag),
                       shrinkflg=bool(pre_shrink_flag),
                       overrideflg=bool(force),
                       rawscflg=bool(use_raw_not_ppv),
                       apcflg=bool(apply_apc_flg),
                       rhodefault=float(set_default_rho),
                       targfnzero=float(target_fraction_of_none_zero),
                       thresh=float(set_lasso_convergence_threshold),
                       idthresh=float(identity_threshold),
                       pseudoc=int(set_pseudocount_value),
                       minseqsep=int(set_minimum_sequence_separation),
                       blockfn=str(rho_parameter_file),
                       maxgapf=float(maximum_fraction_of_gaps),
                       maxthread=int(maxthread))

    if isinstance(psicov, tuple) and psicov[0] == None:
        if psicov[1] == 0:
            return psicov
        elif psicov[1] == 1:
            from ..IO.output import printError
            printError("Out of memory!")
            return None
        elif psicov[1] == 2:
            from ..IO.output import printError
            printError(
                "Not enough sequences or sequence diversity to proceed!")
            printError(
                "Neff ({0}) < MINEFSEQS ({1})".format(psicov[2], msa.shape[1]))
            printError(
                "If you want to force a calculation at your own risk, use force=1.")
            return None

    return psicov


def buildPSICOV(msa,
                approx_lasso_flag=0,
                force=0,
                target_fraction_of_none_zero=0,
                set_minimum_sequence_separation=5,
                maxthread=-1,
                **kwargs):
    """Return PSICOV matrix calculated for *msa*, which may be an
    :class:`.MSA` instance or a 2D Numpy character array.

    You could use `buildPSICOV_expert` for more detail usage.

    Options:
        `approx_lasso_flag`: use approximate LASSO algorithm, default as 0.
        `force`: force the calculation if the sequence is not sufficient,
            default as 0.
        `target_fraction_of_none_zero`: set the target matrix density value 
            (none-zero fraction), should be in range 5e-5 - 1. default as 0
            (not specified)
        `set_minimum_sequence_separation`: set pseudocount value, default as 5.
        `maxthread`: the max number of threads, default as -1(choose max).

    [JDT12] Jones, David T., et al. "PSICOV: precise structural contact
    prediction using sparse inverse covariance estimation on large multiple
    sequence alignments." Bioinformatics 28.2 (2012): 184-190."""

    msa = getMSA(msa)

    pre_shrink_flag = 1
    use_raw_not_ppv = 1
    apply_apc_flg = 1
    set_default_rho = -1
    set_lasso_convergence_threshold = 1e-4
    identity_threshold = -1
    set_pseudocount_value = 1
    rho_parameter_file = ""
    maximum_fraction_of_gaps = 0.9

    if (not 5e-5 < float(target_fraction_of_none_zero) < 1) and (target_fraction_of_none_zero != 0):
        from ..IO.output import printError
        printError("target_fraction_of_none_zero must between 5e-5 and 1.")
        return None
    from numpy import zeros
    from .Ccorrelation_p import msapsicov
    length = msa.shape[1]
    psicov = zeros((length, length), float)
    psicov = msapsicov(msa, psicov,
                       approxflg=bool(approx_lasso_flag),
                       shrinkflg=bool(pre_shrink_flag),
                       overrideflg=bool(force),
                       rawscflg=bool(use_raw_not_ppv),
                       apcflg=bool(apply_apc_flg),
                       rhodefault=float(set_default_rho),
                       targfnzero=float(target_fraction_of_none_zero),
                       thresh=float(set_lasso_convergence_threshold),
                       idthresh=float(identity_threshold),
                       pseudoc=int(set_pseudocount_value),
                       minseqsep=int(set_minimum_sequence_separation),
                       blockfn=str(rho_parameter_file),
                       maxgapf=float(maximum_fraction_of_gaps),
                       maxthread=int(maxthread))

    if isinstance(psicov, tuple) and psicov[0] == None:
        if psicov[1] == 0:
            return psicov
        elif psicov[1] == 1:
            from ..IO.output import printError
            printError("Out of memory!")
            return None
        elif psicov[1] == 2:
            from ..IO.output import printError
            printError(
                "Not enough sequences or sequence diversity to proceed!")
            printError(
                "Neff ({0}) < MINEFSEQS ({1})".format(psicov[2], msa.shape[1]))
            printError(
                "If you want to force a calculation at your own risk, use force=1.")
            return None

    return psicov


def applyPPV(psicov, **kwargs):
    """Apply the PPV test for PSICOV. If you have used `use_raw_not_ppv`=0 in
    `buildPSICOV_expert`, you could skip this step. The PPV model was built in
    [TN12]. The PPV refinement scales the PSICOV scores to 0-0.904 (which is
    the positive predictive values in statistics).

    [TN12] Nugent T, Jones D T. Accurate de novo structure prediction of large
    transmembrane protein domains using fragment-assembly and correlated mutation
    analysis[J]. Proceedings of the National Academy of Sciences, 2012,109(24):
    E1540-E1547."""

    from ..IO.output import printInfo
    printInfo("The given matrix should be the PSICOV matrix with no `use_raw_not_ppv`=0 "
              "option in `buildPSICOV_expert` or any from `buildPSICOV`.")

    from numpy import e, zeros_like

    r = zeros_like(psicov)
    r[psicov != 0] = 0.904 / (1. +
                              16.61 * (e**(psicov[psicov != 0] * (-0.8105))))
    return r
