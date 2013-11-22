'''Some Fasta file read and write functions.
'''

__author__ = 'Wenzhi Mao'

__all__ = ['ReadFasta']


def _Startup():
    '''Get _path__.'''
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


def ReadFasta(filename, title=False, **kwargs):
    '''This a function to read Fasta file.
    Given a filename and the function will return a list of sequences.
    If the Fasta sequences name are requested, add Name=True to get name also.'''
    from numpy import array
    f = open(filename, 'r')
    fastas = f.read()[1:].split('>')
    f.close()
    name = array([fasta.split('\n')[0] for fasta in fastas])
    seq = array([list(''.join(''.join(fasta.split('\n')[1:]).split())) for fasta in fastas])
    if title:
        return name, seq
    else:
        return seq


_Startup()
