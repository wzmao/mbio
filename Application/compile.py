'''It is a module to compile C files.'''


__author__ = 'Wenzhi Mao'
__all__ = []


def _Startup():
    '''Get _path__.'''
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


def make(p, compiler='', option=''):
    '''Compile C files using mpicc.'''
    from os import path
    from os import popen
    abp = path.abspath(p)
    if path.splitext(path.split(abp)[1])[1] == '.c':
        sop = abp[:-2]+'_c.so'
    compiler = compiler if compiler else 'gcc'
    option = option if option else '-shared -fPIC -O3'
    if popen('which {0}'.format(compiler)).read()=='':
        from mbio.IO import error
        print 'adfad'
        error.errorprint('The compler `{0}` is not available on this computer. gcc used instead.'.format(compiler))
        compiler='gcc'
    popen('{0} {1} -o {2} {3}'.format(compiler, option, sop, abp))