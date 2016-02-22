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
           'buildPSICOV',
           'calcMeff', 'calcContactFrac',
           'applyAPC', 'applyBND', 'applyPPV',
           'applyDICOV']


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
    from ..Application.math import invsp

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

    c = invsp(c)

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
    from ..Application.math import eigh
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
    d, u = eigh(mat_th)
    for i in xrange(n):
        if d[i] != 0:
            d[i] = (-1. + (1 + 4 * d[i] * d[i]) ** .5) / 2. / d[i]
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
                              16.61 * (e ** (psicov[psicov != 0] * (-0.8105))))
    return r


def calcContactFrac(n, **kwargs):
    """Return the proportion of contact for a n-residue protein.
    This proportion is utilized in the DICOV as the prior probability of 3D
    contact, P(+), by a regression analysis of a training set of 162 structurally
    known protein sequences.

    [MW15] Mao W, Kaya C, Dutta A, et al. Comparative study of the effectiveness
    and limitations of current methods for detecting sequence coevolution[J].
    Bioinformatics, 2015: btv103."""

    return 8.88332505338 * (n ** -0.825688603083)


def applyDICOV(msa=None, di=None, dca=None, psicov=None, **kwargs):
    """This function use a naïve Bayes classifier to combine DCA and PSICOV as
    described in [MW15]. The calculation is based on DI and PSICOV. You could
    provide these two matrices or the MSA. The PSICOV matrix should be the PPV
    scaled. You could use `buildPSICOV_expert` and use `use_raw_not_ppv`=0 option
    or use `applyPPV` to build it from the pre-PPV PSICOV matrix.

    [MW15] Mao W, Kaya C, Dutta A, et al. Comparative study of the effectiveness
    and limitations of current methods for detecting sequence coevolution[J].
    Bioinformatics, 2015: btv103."""

    from numpy import zeros_like, log, mgrid, fromfile, array, cov, diag, e, pi, trunc
    from numpy.linalg.linalg import inv
    from ..IO.output import printError
    from ..Constant import getconstantfunc

    if ((di == None and dca == None) or (psicov == None)) and msa == None:
        printError("DI and PSICOV matrices or MSA should be provided.")
        return None
    if di != None and dca != None and not (di == dca).all():
        printError(
            "DI and DCA matrices are not the same, check it or just use one.")
        return None

    d = di if di != None else dca if dca != None else None
    p = psicov if psicov != None else None

    if d != None:
        if d.ndim != 2:
            printError("The dimension of DI matrix is wrong.")
            return None
        elif d.shape[0] != d.shape[1]:
            printError("DI matrix is not square.")
            return None
        elif p == None and getMSA(msa).shape[1] != d.shape[0]:
            printError("DI matrix does not fit the MSA.")
            return None
    if p != None:
        if p.ndim != 2:
            printError("The dimension of PSICOV matrix is wrong.")
            return None
        elif p.shape[0] != p.shape[1]:
            printError("PSICOV matrix is not square.")
            return None
        elif d == None and getMSA(msa).shape[1] != p.shape[0]:
            printError("PSICOV matrix does not fit the MSA.")
            return None
    if p != None and d != None:
        if p.shape[0] != d.shape[0]:
            printError("DI and PSICOV matrices have different sizes.")
            return None

    d = buildDI(msa) if d == None else d
    p = buildPSICOV_expert(msa, use_raw_not_ppv=0) if p == None else p

    n = di.shape[0]

    dicov = zeros_like(d)
    pplus = getconstantfunc('pplus')
    pminus = getconstantfunc('pminus')
    X, Y = mgrid[-5:0.05:0.05, -3:0.05:0.05]
    pplus.resize(X.shape)
    pminus.resize(X.shape)

    psigplus = getconstantfunc('psigplus')
    psigminus = getconstantfunc('psigminus')
    XX = mgrid[-5:0.05:0.05]
    psigplus.resize(XX.shape)
    psigminus.resize(XX.shape)

    prate = calcContactFrac(d.shape[0])
    qrate = 1.0 - prate

    pplus = pplus * prate
    psigplus = psigplus * prate
    pminus = pminus * qrate
    psigminus = psigminus * qrate

    cdouble = zeros_like(pplus)

    pos = []
    val = []
    for i in xrange(pplus.shape[0]):
        for j in xrange(pplus.shape[1]):
            if pplus[i][j] <= 1e-3 or pminus[i][j] <= 1e-3:
                cdouble[i][j] = -1
            else:
                cdouble[i][j] = pplus[i][j] / (pplus[i][j] + pminus[i][j])
                pos.append([X[i][j], Y[i][j]])
                val.append(cdouble[i][j])
    pos = array(pos).T
    val = array(val)
    nn = val.shape[0]
    s = cov(pos)
    invs = inv(s)
    h = (nn ** (-1.0 / 6.0)) * (((2.0 ** (-1.0)) * (diag(s).sum())) ** .5)
    para1 = (-.5) * (h ** -2.)
    for i in xrange(pplus.shape[0]):
        for j in xrange(pplus.shape[1]):
            if cdouble[i][j] == -1:
                temp = array([X[i][j], Y[i][j]]).reshape((2, 1))
                temp = e ** (para1 *
                            ((invs.dot((pos - temp)) * (pos - temp)).sum(0)))
                cdouble[i][j] = (temp.dot(val).sum()) / temp.sum()

    csingle = zeros_like(psigplus)
    pos = []
    val = []
    for i in xrange(psigplus.shape[0]):
        if psigplus[i] <= 1e-8 or psigminus[i] <= 1e-8:
            csingle[i] = -1
        else:
            csingle[i] = psigplus[i] / (psigplus[i] + psigminus[i])
            pos.append(XX[i])
            val.append(csingle[i])
    pos = array(pos).T
    val = array(val)
    nn = val.shape[0]
    s = pos.std()
    h = ((4.0 / 3) ** (.2)) * s * (nn ** -.2)
    para1 = 1. / (nn * h * ((2 * pi) ** .5))
    for i in xrange(psigplus.shape[0]):
        if csingle[i] == -1:
            temp = para1 * (e ** (-.5 * (((pos - XX[i]) / h) ** 2)))
            csingle[i] = (temp.dot(val).sum()) / temp.sum()

    for i in xrange(n):
        for j in xrange(i + 1, n):
            if psicov[i][j] == 0:
                ldi = log(di[i][j]) / log(10)
                temp = int(trunc((ldi + 5) / .05))
                dicov[i][j] = dicov[j, i] = csingle[
                    temp] + (ldi - XX[temp]) / .05 * (csingle[temp + 1] - csingle[temp])
            else:
                ldi = log(di[i][j]) / log(10)
                lps = log(psicov[i][j]) / log(10)
                temp1 = int(trunc((ldi + 5) / .05))
                temp2 = int(trunc((lps + 3) / .05))
                val = array([[cdouble[temp1, temp2], cdouble[
                            temp1, temp2 + 1]], [cdouble[temp1 + 1, temp2], cdouble[temp1 + 1, temp2 + 1]]])
                dicov[i, j] = dicov[j, i] = array([X[temp1 + 1][temp2] - ldi, ldi - X[temp1][temp2]]).dot(
                    val).dot(array([Y[temp1][temp2 + 1] - lps, lps - Y[temp1][temp2]])) * 400

    return dicov

# def buildplmDCA(msa,refine=0,unique=0,weighting=0.2,maxthread=-1,**kwargs):

#     msa = getMSA(msa)

#     from numpy import empty,ones
#     from .Ccorrelation import msaplmdca

#     if refine:
#         mark = msa[0]!='-'
#         for i in xrange(len(mark)):
#             if mark[i] and msa[0,i].isupper():
#                 mark[i]=True
#             else:
#                 mark[i]=False
#         msa=msa[:,mark]
#     if unique:
#         mark=ones(msa.shape[0],dtype=bool)
#         for i in xrange(msa.shape[0]):
#             if not mark[i]:
#                 continue
#             for j in xrange(i+1,msa.shape[0]):
#                 if (msa[i]==msa[j]).all():
#                     mark[j]=False
#         msa=msa[mark]

#     length = msa.shape[1]
#     plmdca = empty((length, length), float)
#     plmdca = msaplmdca(msa, plmdca, weighting=weighting)

#     return plmdca
