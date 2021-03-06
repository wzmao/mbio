# -*- coding: utf-8 -*-
"""MSA class to store MSA sequence.
"""

__author__ = 'Wenzhi Mao'
__all__ = ['MSA']


class MSA(object):

    """A MSA class to store sequences indexed by labels."""

    def __init__(self, msa, labels=None, **kwargs):
        from numpy import dtype, array
        if type(msa) == list:
            msa = [list(i) for i in msa]
            msa = array(msa, dtype='|S1')
        if not set(['ndim', 'shape']).issubset(set(dir(msa))):
            raise TypeError('MSA must be numpy array or list')

        if msa.ndim != 2:
            raise ValueError('Must be 2-dimension.')
        (numseq, numres) = msa.shape

        if labels == None:
            labels = array(['None'] * numseq, dtype='|S4')
        if type(labels) == list:
            labels = array(labels)
        if len(labels) != numseq:
            raise ValueError(
                'Label number wrong. {0} Sequences and {1} labels'.format(numseq, len(labels)))

        self.seq = msa
        self.label = labels
        self._dict = {}
        for i in xrange(numseq):
            if not labels[i] in self._dict.keys():
                self._dict[labels[i]] = [i]
            else:
                self._dict[labels[i]] = self._dict[labels[i]] + [i]
        self.numseq = numseq
        self.numres = numres

    def __repr__(self):

        return "MSA = {0} sequences with {1} residues.".format(self.numseq, self.numres)

    def __getitem__(self, index):
        from numpy import array
        if type(index) == int:
            return self.seq[index]
        if type(index) == str:
            if index in self._dict.keys():
                return array([self.seq[i] for i in self._dict[index]])
        raise ValueError('Cannot index.')

    def __iter__(self):

        for i in xrange(self.numseq):
            yield self.seq[i]
