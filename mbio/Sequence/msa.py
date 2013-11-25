'''MSA
'''

__author__ = 'Wenzhi Mao'
__all__ = ['MSA', 'ReadFasta']

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
        self.numseq = numseq
        self.numres = numres

    def __repr__(self):

        return "MSA = {0} sequences with {1} residues.".format(self.numseq, self.numres)

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


def ReadFasta(filename, **kwargs):
    '''This a function to read Fasta file.
    Given a filename and the function will return a list of sequences.'''
    from numpy import array
    from os.path import isfile
    from .c_fasta import readFasta
    if not isfile(filename):
        raise ValueError("File not found.")
    msa = readFasta(filename)
    msa = MSA(msa[1], labels=msa[0])
    return msa

_Startup()
