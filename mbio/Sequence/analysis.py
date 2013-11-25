'''Some MSA and correlation matrix analysis functions.
'''

__author__ = 'Wenzhi Mao'
__all__ = ['CombineMSA']

from numpy import dtype, array


def _Startup():
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()

def CombineMSA(msa1, msa2, spi=1, prot=0, **kwargs):
	from numpy import concatenate
	# try:
	# 	seq=concatenate(msa1.seq,msa2.seq,ax)
	return 0


_Startup()
