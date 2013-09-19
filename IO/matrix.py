'''Some binary matrix read functions.
'''

__author__ = 'Wenzhi Mao'

__all__ = ['ReadMatrix']


def _Startup():
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


def ReadMatrix(filename, dtype, l=None):
    import struct
    f = open(filename, 'r')
    tempmatrix = f.read()
    f.close()
    if l:
        if len(tempmatrix) != {'d': 8, 'i': 4}[dtype]*(l**2):
            print '* OMES file size wrong.'
            return None
    else:
        l = int((len(tempmatrix)/{'d': 8, 'i': 4}[dtype])**.5)
    m = []
    for i in range(l):
        m.append(list(struct.unpack(str(l)+dtype, tempmatrix[i*l*{
                 'd': 8, 'i': 4}[dtype]:(i*l+l)*{'d': 8, 'i': 4}[dtype]])))
    return m


_Startup()
