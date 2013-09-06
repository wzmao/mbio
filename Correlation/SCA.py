'''SCA(statistical coupling analysis) is a correlation algrithm.
'''
__author__ = 'Wenzhi Mao'
__all__ = ['CalcSCA', 'CalcSCA1']


def _Startup():
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()
    _parse()


def _parse():
    pass


def CalcSCA(sequence):
    from math import fsum as sum
    from math import log

    def mean(x):
        return 1.0*sum(x)/len(x)
    reslist = [i for i in 'ACDEFGHIKLMNPQRSTVWY']
    q = [.073, .025, .050, .061, .042, .072, .023, .053, .064, .089,
         .023, .043, .052, .040, .052, .073, .056, .063, .013, .033]
    p = []
    phi = []
    m = len(sequence)
    l = len(sequence[0])
    for i in range(l):
        temp = [sequence[j][i] for j in range(m)]
        f = [temp.count(j)*1.0/m for j in reslist]
        p.append([])
        phi.append([])
        sumphif = 0
        for j in range(len(q)):
            if 0 < f[j] < 1:
                phi[-1].append(abs(log(f[j]*(1-q[j])/(1-f[j])/q[j])))
                p[-1].append(f[j]*phi[-1][-1])
            else:
                phi[-1].append(0)
                p[-1].append(0)
            sumphif += p[-1][-1]**2
        sumphif = sumphif**.5
        if sumphif == 0.0:
            for j in range(len(p[-1])):
                p[-1][j] = 0.0
        else:
            for j in range(len(p[-1])):
                p[-1][j] /= sumphif
    xp = []
    for i in range(m):
        xp.append([])
        for j in range(l):
            if sequence[i][j] in reslist:
                xp[-1].append(p[j][reslist.index(sequence[i][
                              j])]*phi[j][reslist.index(sequence[i][j])])
            else:
                xp[-1].append(0)
    meanxp = [mean([xp[j][i] for j in range(m)]) for i in range(l)]
    c = [[0 for j in range(l)]for i in range(l)]
    for i in range(l):
        for j in range(i, l):
            c[i][j] = abs(sum([xp[k][i]*xp[k][
                          j] for k in range(m)])*1.0/m-meanxp[i]*meanxp[j])
            c[j][i] = c[i][j]
    return c,xp

def CalcSCA1(sequence):
    from math import fsum as sum
    from math import log

    def mean(x):
        return 1.0*sum(x)/len(x)
    reslist = [i for i in 'ACDEFGHIKLMNPQRSTVWY']
    q = [.073, .025, .050, .061, .042, .072, .023, .053, .064, .089,
         .023, .043, .052, .040, .052, .073, .056, .063, .013, .033]
    p = []
    phi = []
    m = len(sequence)
    l = len(sequence[0])
    for i in range(l):
        temp = [sequence[j][i] for j in range(m)]
        f = [temp.count(j)*1.0/m for j in reslist]
        p.append([])
        phi.append([])
        sumphif = 0
        for j in range(len(q)):
            if 0 < f[j] < 1:
                phi[-1].append(abs(log(f[j]*(1-q[j])/(1-f[j])/q[j])))
                p[-1].append(f[j]*phi[-1][-1])
            else:
                phi[-1].append(0)
                p[-1].append(0)
            sumphif += p[-1][-1]**2
        sumphif = sumphif**.5
        if sumphif == 0.0:
            for j in range(len(p[-1])):
                p[-1][j] = 0.0
        else:
            for j in range(len(p[-1])):
                p[-1][j] /= sumphif
    xp = []
    for i in range(m):
        xp.append([])
        for j in range(l):
            if sequence[i][j] in reslist:
                xp[-1].append(p[j][reslist.index(sequence[i][
                              j])]*phi[j][reslist.index(sequence[i][j])])
            else:
                xp[-1].append(0)
    meanxp = [mean([xp[j][i] for j in range(m)]) for i in range(l)]
    c = [[0 for j in range(l)]for i in range(l)]
    for i in range(l):
        for j in range(i, l):
            c[i][j] = abs(sum([xp[k][i]*xp[k][
                          j] for k in range(m)])*1.0/m-meanxp[i]*meanxp[j])
            c[j][i] = c[i][j]
    return c,xp

_Startup()