"""This module contains some MSA asistance functions.
"""

__author__ = 'Wenzhi Mao'
__all__ = ['combineMSA']


def combineMSA(msas, spi=1, prot=0, **kwargs):
    """Combine 2 or more MSAs by the labels."""

    from numpy import concatenate, dtype, array
    from numpy.core.shape_base import vstack
    from .msa import MSA
    try:
        n = len(msas)
        labels = array([msas[i].label for i in range(n)])
        numress = array([msas[i].numres for i in range(n)])
    except:
        raise ValueError("They are not MSAs list or array.")
    if ((spi) and (prot)):
        dicts = [{} for i in range(n)]
        for i in range(n):
            for j in range(len(labels[i])):
                dicts[i][labels[i][j].split('/')[0]] = j
        label = set(dicts[0].keys())
        for i in range(1, n):
            label = label.intersection(set(dicts[i].keys()))
        label = list(label)
        # label.sort()
        if len(label) != 0:
            seq = array([concatenate(([msas[i][dicts[i][label[
                        j]]] for i in range(n)])) for j in range(len(label))])
        else:
            return None
        for i in range(len(label)):
            temp = ''
            for j in range(n):
                temp += '+' + msas[j].label[dicts[j][label[i]]].split('/')[1]
            label[i] += '/' + temp.strip('+')
        label = array(label)
        newMSA = MSA(seq, label)
    elif ((spi) and (not prot)):
        freqs = [[(msas[i].seq[:, j] == '-').sum() * 1. / msas[
                  i].numseq for j in range(msas[i].numres)]for i in range(n)]
        score = [[0 for j in range(msas[i].numseq)]for i in range(n)]
        for i in range(n):
            for j in range(msas[i].numseq):
                count = 0.
                for k in range(len(msas[i][j])):
                    if msas[i][j][k] != '-':
                        score[i][j] += freqs[i][k]
                        count += 1.
                score[i][j] /= count
        dicts = [{} for i in range(n)]
        for i in range(n):
            for j in range(len(labels[i])):
                if labels[i][j].split('/')[0].split('_')[1] in dicts[i].keys() and score[i][dicts[i][labels[i][j].split('/')[0].split('_')[1]]] < score[i][j]:
                    dicts[i][labels[i][j].split('/')[0].split('_')[1]] = j
                else:
                    dicts[i][labels[i][j].split('/')[0].split('_')[1]] = j
        label = set(dicts[0].keys())
        for i in range(1, n):
            label = label.intersection(set(dicts[i].keys()))
        label = list(label)
        label.sort()
        if len(label) != 0:
            seq = array([concatenate(([msas[i][dicts[i][label[
                        j]]] for i in range(n)])) for j in range(len(label))])
        else:
            return None
        for i in range(len(label)):
            temp = ''
            for j in range(n):
                temp += '---' + msas[j].label[dicts[j][label[i]]]
            label[i] = temp.strip('---')
        label = array(label)
        newMSA = MSA(seq, label)
    return newMSA
