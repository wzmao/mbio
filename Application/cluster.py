__author__ = 'Wenzhi Mao'
__all__ = []


def _Startup():
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


def Iscluster():
    '''The function will return a True or False for the cluster availability.'''
    from os import popen
    if popen('which qsub').read() == '':
        return False
    return True


def Getclusterlist(available=False):
    '''The function will return a list or None for a list of cluster.
    available option could give assigned to be True for available cluster.'''
    from mbio import _hostname
    if _hostname == 'n000':
        from os import popen
        clist = popen('qstat -f').readlines()
        clist = [i for i in clist if i.startswith('parallel') and
                 i.split()[2].count('/') == 2 and i.split()[3].count('.') != 0]
        clist = [[i.split()[0].split('@')[1].split('.')[0], i.split()[2]]
                         for i in clist if i.split()[2].split('/')[-1] != 0]
        if not available:
            return clist
        else:
            clist = [i for i in clist if
                     i[1].split('/')[-1] != i[1].split('/')[-2]]
            return clist
    else:
        return None


_Startup()
