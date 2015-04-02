'''Some error and information print.'''

__author__ = 'Wenzhi Mao'

__all__ = ['errorprint', 'infoprint']


def errorprint(x):
    '''It is a print red error information function.'''
    print '\x1b[1;31m* {0}\x1b[0;29m'.format(str(x))


def infoprint(x):
    '''It is a print normal information function.'''
    print '\x1b[0;29m* {0}\x1b[0;29m'.format(str(x))
