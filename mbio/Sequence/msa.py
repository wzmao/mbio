'''MSA
'''

__author__ = 'Wenzhi Mao'
__all__ = ['MSA']

from numpy import dtype, array


def _Startup():
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


class MSA(object):

    def __init__(self, msa, labels=None, **kwargs):
        if type(msa) == list:
            msa = [list(i) for i in msa]
            msa = array(msa, dtype='|S1')
        if not set(['ndim', 'shape']).issubset(set(dir(msa))):
            raise TypeError('MSA must be numpy array or list')

        if msa.ndim != 2:
            raise ValueError('Must be 2-dimension.')
        numseq, numres = msa.shape

        if labels == None:
            labels = array(['None'] * numseq, dtype='|S4')
        if type(labels) == list:
            labels = array(labels)
        if len(labels) != numseq:
            raise ValueError(
                'Label number wrong. {} Sequences and {} labels'.format(numseq, len(labels)))

        self.seq = msa
        self.label = labels
        self.numseq = numseq
        self.numres = numres

    def __repr__(self):

        return "MSA = {0} sequences with {0} residues.".format(self.numseq, self.numres)

    def __getitem__(self, index):
        if type(index) == int:
            return self.seq[index]
        if type(index) == str:
            if str in self.label:
                return self.seq[self.label.index(index)]
        raise ValueError('Cannot index.')

    def __iter__(self):

        for i in range(self.numseq):
            yield self.seq[i]

_Startup()
