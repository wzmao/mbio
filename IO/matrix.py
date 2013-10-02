'''Some binary matrix read functions.
'''

__author__ = 'Wenzhi Mao'

__all__ = ['ReadMatrix', 'WriteMatrix']


def _Startup():
    '''Get _path__.'''
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


def ReadMatrix(filename, dtype=None, l=None):
    '''Read Matrix if possible.'''
    if (isinstance(filename, (tuple, list))):
        if len(filename) == 1:
            filename = filename[0]
        elif len(filename) > 1:
            filename = list(filename)[:2]
        else:
            filename = None
    elif type(filename) == str:
        filename = filename
    else:
        filename = None
    if not filename:
        return None
    elif type(filename) == list:
        mlist = [parsesinglematrix(i) for i in filename]
        if None in mlist:
            from mbio.IO import error
            error.errorprint('There is non-square matrix.')
            return None
        llist = [mlist[0][2], mlist[1][2]]
        if len(set(llist)) != 1 or (l and len(set(llist+[l])) != 1):
            from mbio.IO import error
            error.errorprint('Sizes of matrix are different.')
            return None
        else:
            return [mlist[0][0], mlist[1][0]]
    else:
        m = parsesinglematrix(filename)
        if m == None:
            from mbio.IO import error
            error.errorprint('It is not a square matrix.')
            return None
        if l and l != m[2]:
            from mbio.IO import error
            error.errorprint('Sizes of matrix are different.')
            return None
        if dtype and dtype != m[1]:
            from mbio.IO import error
            error.errorprint('Data type are different.')
            return None
        else:
            return m[0]


def parsesinglematrix(filename):
    import struct
    from mbio.Application import math as mmath
    f = open(filename, 'r')
    tempmatrix = f.read()
    f.close()
    if mmath.issquare(len(tempmatrix)/4):
        dtype = 'i'
    elif mmath.issquare(len(tempmatrix)/8):
        dtype = 'd'
    else:
        return None
    l = int((len(tempmatrix)/{'d': 8, 'i': 4}[dtype])**.5)
    m = []
    for i in range(l):
        m.append(list(struct.unpack(str(l)+dtype, tempmatrix[i*l*{
                 'd': 8, 'i': 4}[dtype]:(i*l+l)*{'d': 8, 'i': 4}[dtype]])))
    if dtype == 'i':
        mm = [[0.0 for j in range(l)]for i in range(l)]
        for i in range(l):
            for j in range(l):
                mm[i][j] = m[i][j]+m[j][i]
        m = mm
    return [m, dtype, l]


def WriteMatrix(m, filename, dtype=''):
    l = list(set([len(m)]+[len(i) for i in m]))
    if len(l) != 1:
        from mbio.IO import error
        error.errorprint('The matrix is not square.')
        return None
    l = l[0]
    if not dtype or not dtype in ['i', 'd']:
        dtype = {int: 'i', float: 'd'}[type(m[0][0])]
    import struct
    f = open(filename, 'w')
    st = ''
    for ii in m:
        for jj in ii:
            st += struct.pack(dtype, jj)
    f.write(st)
    f.close()


_Startup()
