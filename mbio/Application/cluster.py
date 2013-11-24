'''It is a module could get some information about current computer.
Sun Grid Engine could be monitored by this module.'''

__author__ = 'Wenzhi Mao'
__all__ = []


def _Startup():
    '''Get _path__.'''
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


def Iscluster():
    '''The function will return a True or False for the cluster availability.'''
    from os import popen
    if popen('which qstat').read() == '':
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
                 for i in clist if i.split()[2].split('/')[-1] != '0']
        if not available:
            return clist
        else:
            clist = [i for i in clist if
                     i[1].split('/')[-1] != i[1].split('/')[-2]]
            return clist
    else:
        return None


def Getname(minc=0, maxc=64):
    '''The function return the cluster name and the available nodes number.
    min and max of nodes number could be assigned.'''
    clist = Getclusterlist(True)
    clist = [[int(i[1].split('/')[-1]) - int(i[1].split('/')[-2])] + i for i in clist
             if minc <= int(i[1].split('/')[-1]) - int(i[1].split('/')[-2]) <= maxc]
    clist.sort(reverse=True)
    if len(clist) == 0:
        return []
    else:
        if clist[0][2].split('/')[-1] != '4':
            result = [{'64': '@bahar64',
                       '32': '@core32',
                       '24': '@core24',
                       '12': '@bahar12',
                       '8': '@core8',
                       }[clist[0][2].split('/')[-1]]]
        else:
            if 92 <= int(clist[0][1][1:]) <= 102:
                result = '@bahar4'
            else:
                result = '@core4'
        return [result, str(clist[0][0])]

_Startup()
