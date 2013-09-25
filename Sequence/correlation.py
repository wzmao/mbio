'''Different algorithm could identify protein position correlation.
We have MI, MIp, OMES and SCA now.
'''

__author__ = 'Wenzhi Mao'
__all__ = ['CalcMI', 'CalcMIp', 'CalcOMES', 'CalcSCA', 'CalcMy']


def _Startup():
    '''Get _path__.'''
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


def CalcMI(sequences):
    '''It is a function to calculate the MI matrix based on language C.
    Given the sequences in a list with no format.
    '''
    if not '_c_CalcMI' in globals().keys():
        import ctypes as ct
        from os import path
        M = ct.CDLL(path.join(_path__, 'correlation_c.so'))
        global _c_CalcMI
        _c_CalcMI = M.calcMI
        _c_CalcMI.argtypes = [ct.POINTER(ct.c_char), ct.c_int, ct.c_int]
        _c_CalcMI.restype = ct.POINTER(ct.c_double)
        _c_CalcMI.__doc__ = '''It is a function to calulate the MI matrix in C.
        Give 4 variable.    `M, n, l`.
            M is the sequence array with length n*l.(char array)
            n is the sequence number.(int)
            l is the sequence length.(int)
        Return 1 variable.  `mi`
            mi is an array with length l*l to return the result.(double array)'''
    import ctypes as ct
    allsequence = ''.join(sequences)
    m = (ct.c_char * len(allsequence))()
    for i in range(len(allsequence)):
        m[i] = allsequence[i]
    l = len(sequences[0])
    result = _c_CalcMI(m, len(sequences), l)
    mi = []
    for i in range(l**2):
        if i % l == 0:
            mi.append([])
        mi[-1].append(result[i])
    return mi


def CalcMIp(sequences):
    '''It is a function to calculate the MIp matrix based on language C.
    Given the sequences in a list with no format.
    '''
    # TODO it will return nan instead of 0.0 for one sequence.
    if not '_c_CalcMIp' in globals().keys():
        import ctypes as ct
        from os import path
        M = ct.CDLL(path.join(_path__, 'correlation_c.so'))
        global _c_CalcMIp
        _c_CalcMIp = M.calcMIp
        _c_CalcMIp.argtypes = [ct.POINTER(ct.c_char), ct.c_int, ct.c_int]
        _c_CalcMIp.restype = ct.POINTER(ct.c_double)
        _c_CalcMIp.__doc__ = '''It is a function to calulate the MIp matrix in C.
        Give 4 variable.    `M, n, l`.
            M is the sequence array with length n*l.(char array)
            n is the sequence number.(int)
            l is the sequence length.(int)
        Return 1 variable.  `mip`
            mip is an array with length l*l to return the result.(double array)'''
    import ctypes as ct
    allsequence = ''.join(sequences)
    m = (ct.c_char * len(allsequence))()
    for i in range(len(allsequence)):
        m[i] = allsequence[i]
    l = len(sequences[0])
    result = _c_CalcMIp(m, len(sequences), l)
    mip = []
    for i in range(l**2):
        if i % l == 0:
            mip.append([])
        mip[-1].append(result[i])
    return mip


def CalcOMES(sequences):
    '''It is a function to calculate the OMES matrix based on language C.
    Given the sequences in a list with no format.
    '''
    if not '_c_CalcOMES' in globals().keys():
        import ctypes as ct
        from os import path
        M = ct.CDLL(path.join(_path__, 'correlation_c.so'))
        global _c_CalcOMES
        _c_CalcOMES = M.calcOMES
        _c_CalcOMES.argtypes = [ct.POINTER(ct.c_char), ct.c_int, ct.c_int]
        _c_CalcOMES.restype = ct.POINTER(ct.c_double)
        _c_CalcOMES.__doc__ = '''It is a function to calulate the OMES matrix in C.
        Give 4 variable.    `M, n, l`.
            M is the sequence array with length n*l.(char array)
            n is the sequence number.(int)
            l is the sequence length.(int)
        Return 1 variable.  `omes`
            omes is an array with length l*l to return the result.(double array)'''
    import ctypes as ct
    allsequence = ''.join(sequences)
    m = (ct.c_char * len(allsequence))()
    for i in range(len(allsequence)):
        m[i] = allsequence[i]
    l = len(sequences[0])
    result = _c_CalcOMES(m, len(sequences), l)
    omes = []
    for i in range(l**2):
        if i % l == 0:
            omes.append([])
        omes[-1].append(result[i])
    return omes


def CalcSCA(sequences):
    '''It is a function to calculate the SCA matrix based on language C.
    Given the sequences in a list with no format.
    '''
    if not '_c_CalcSCA' in globals().keys():
        import ctypes as ct
        from os import path
        M = ct.CDLL(path.join(_path__, 'correlation_c.so'))
        global _c_CalcSCA
        _c_CalcSCA = M.calcSCA
        _c_CalcSCA.argtypes = [ct.POINTER(ct.c_char), ct.c_int, ct.c_int]
        _c_CalcSCA.restype = ct.POINTER(ct.c_double)
        _c_CalcSCA.__doc__ = '''It is a function to calulate the SCA matrix in C.
        Give 4 variable.    `M, n, l`.
            M is the sequence array with length n*l.(char array)
            n is the sequence number.(int)
            l is the sequence length.(int)
        Return 1 variable.  `sca`
            sca is an array with length l*l to return the result.(double array)'''
    import ctypes as ct
    allsequence = ''.join(sequences)
    m = (ct.c_char * len(allsequence))()
    for i in range(len(allsequence)):
        m[i] = allsequence[i]
    l = len(sequences[0])
    sca = []
    result = _c_CalcSCA(m, len(sequences), l)
    for i in range(l**2):
        if i % l == 0:
            sca.append([])
        sca[-1].append(result[i])
    return sca


'''
def CalcSCA1(sequence):
    from math import fsum as sum
    from math import log

    def mean(x):
        return 1.0*sum(x)/len(x)
    reslist = [i for i in 'ACDEFGHIKLMNPQRSTVWY']
    q = [.073, .025, .050, .061, .042, .072, .023, .053, .064, .089,
         .023, .043, .052, .040, .052, .073, .056, .063, .013, .033]
    p = []
    phi = []
    m = len(sequence)
    l = len(sequence[0])
    for i in range(l):
        temp = [sequence[j][i] for j in range(m)]
        f = [temp.count(j)*1.0/m for j in reslist]
        p.append([])
        phi.append([])
        sumphif = 0
        for j in range(len(q)):
            if 0 < f[j] < 1:
                phi[-1].append(abs(log(f[j]*(1-q[j])/(1-f[j])/q[j])))
                p[-1].append(f[j]*phi[-1][-1])
            else:
                phi[-1].append(0)
                p[-1].append(0)
            sumphif += p[-1][-1]**2
        sumphif = sumphif**.5
        if sumphif == 0.0:
            for j in range(len(p[-1])):
                p[-1][j] = 0.0
        else:
            for j in range(len(p[-1])):
                p[-1][j] /= sumphif
    xp = []
    for i in range(m):
        xp.append([])
        for j in range(l):
            if sequence[i][j] in reslist:
                xp[-1].append(p[j][reslist.index(sequence[i][
                              j])]*phi[j][reslist.index(sequence[i][j])])
            else:
                xp[-1].append(0)
    meanxp = [mean([xp[j][i] for j in range(m)]) for i in range(l)]
    c = [[0 for j in range(l)]for i in range(l)]
    for i in range(l):
        for j in range(i, l):
            c[i][j] = abs(sum([xp[k][i]*xp[k][
                          j] for k in range(m)])*1.0/m-meanxp[i]*meanxp[j])
            c[j][i] = c[i][j]
    return c
'''


def CalcMy(sequences):
    if not '_c_CalcMy' in globals().keys():
        import ctypes as ct
        from os import path
        M = ct.CDLL(path.join(_path__, 'correlation_c.so'))
        global _c_CalcMy
        _c_CalcMy = M.calcMy
        _c_CalcMy.argtypes = [ct.POINTER(ct.c_char), ct.c_int, ct.c_int]
        _c_CalcMy.restype = ct.POINTER(ct.c_double)
        _c_CalcMy.__doc__ = '''It is a function to calulate the SCA matrix in C.
        Give 4 variable.    `M, n, l`.
            M is the sequence array with length n*l.(char array)
            n is the sequence number.(int)
            l is the sequence length.(int)
        Return 1 variable.  `sca`
            sca is an array with length l*l to return the result.(double array)'''
    import ctypes as ct
    allsequence = ''.join(sequences)
    m = (ct.c_char * len(allsequence))()
    for i in range(len(allsequence)):
        m[i] = allsequence[i]
    l = len(sequences[0])
    result = _c_CalcMy(m, len(sequences), l)
    # my = []
    # for i in range(l**2):
    #     if i % l == 0:
    #         sca.append([])
    #     sca[-1].append(result[i])
    return result


_Startup()
