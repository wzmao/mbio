'''MI is a simple protein correlation algrithm.
We use MI to calculate the correlation for 20 kinds of amino acid and gaps.(21)
It also provide a shuffling function to calculate the shuffled P-value.
The shuffled P-value calculation is performced by C and MPI.
'''
__author__ = 'Wenzhi Mao'
__all__ = ['CalcMI', 'CalcMIp', 'CalcOMES']


def _Startup():
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
        M = ct.CDLL(path.join(_path__,'mi_c.so'))
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
    if not '_c_CalcMIp' in globals().keys():
        import ctypes as ct
        from os import path
        M = ct.CDLL(_path__+'/mi_c.so')
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
        M = ct.CDLL(_path__+'/omes_c.so')
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


_Startup()
