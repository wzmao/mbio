__author__ = 'Wenzhi Mao'

__all__ = ['QuickSort']


def _Startup():
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


def QuickSort(l):
    '''This is a quick sort function based on language C.
    Given a list of integer or float numbers, the function will return a sorted list.'''
    if not '_c_qs' in globals().keys():
        import ctypes as ct
        from os import path
        M = ct.CDLL(path.join(_path__, 'sort_c.so'))
        global _c_qs
        _c_qs = M.qs
        _c_qs.argtypes = [ct.POINTER(ct.c_double), ct.c_int, ct.c_int]
        _c_qs.restype = ct.c_void_p
        _c_qs.__doc__ = '''It is a function to quick sort an array in C double format.
        Give 3 variable.    `M, l, h`.
            M is the array.(double array)
            l is 0 or the lower cutoff.(int)
            h is len(M)-1 or the up cutoff.(int)'''
    import ctypes as ct
    carray = (ct.c_double * len(l))()
    for i in range(len(l)):
        carray[i] = l[i]
    _c_qs(carray, 0, len(l)-1)
    result = [i for i in carray]
    return result

_Startup()