'''Get job organized in Sun Grid Engine.
It also contains some allover functions for job organization.
'''

__author__ = 'Wenzhi Mao'
__all__ = []


def _Startup():
    '''Get _path__.'''
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


def AskJobNumber():
    '''Get a job number for new job.'''
    import os
    from os import path
    if not os.path.exists(path.join(_path__, '..', '.Cache')):
        os.mkdir(path.join(_path__, '..', '.Cache'))
    plist = os.listdir(path.join(_path__, '..', '.Cache'))
    i = 1
    while os.path.exists(path.join(_path__, '..', '.Cache', str(i)+'.c')):
        i += 1
    return str(i)


def MkdirResult():
    '''Make .Result folder.'''
    from os import path, mkdir
    if not path.exists(path.join(_path__, '..', '.Result')):
        mkdir(path.join(_path__, '..', '.Result'))


def Writejob(jobnumber, scripts):
    '''Write the scripts into a C file in .Cache folder.'''
    from os import path
    f = open(path.join(_path__, '..', '.Cache', jobnumber+'.c'), 'w')
    f.write(scripts)
    f.close()


def Getpronum():
    '''Get process number for the current computer.
    Used for local run jobs.'''
    from os import popen
    return int(popen('nproc').read())


def SubShufflejob(jobnumber, cluster, core):
    '''Submit shuffle jobs.(ShuffleMI, ShuffleMIp, ShuffleOMES)
    cluster could be assignd as True or False. Otherwise the job will run locally.
    core could be a integer for min core number, or a tuple to set a min and max.'''
    from os import popen, path
    from mbio.Application import cluster
    if cluster and cluster.Iscluster():
        if isinstance(core, (list, tuple)) and len(core) != 0:
            cmin = core[0]
            if len(core) > 1:
                cmax = core[1]
            name = Getname(cmin, cmax)
        elif isinstance(core, (int, float)):
            name = Getname(int(core))
        else:
            name = Getname()
        f = open(path.join(_path__, '..', 'Scripts', 'shuffle-qsub.job'), 'r')
        qsub = f.read()
        f.close()
        qsub = [qsub.replace('core_number', name[1]).replace
                ('clustername', name[0]).replace
                ('output_name', jobnumber+'.output').replace
                ('error_name', jobnumber+'.err').replace
                ('c_file_name', jobnumber+'.c').replace
                ('out_file_name', jobnumber+'.out').replace
                ('job_name', 'job_'+jobnumber)][0]
        f = open(path.join(_path__, '..', '.Cache', jobnumber+'.job'), 'w')
        f.write(qsub)
        f.close()
        workid = popen(
            'cd '+path.join(_path__, '..', '.Cache')+';qsub '+jobnumber+'.job').read()
        workid = workid.split()[2]
        import time
        while popen('qstat -j '+str(workid)+' 2>&1').read().startswith('===='):
            time.sleep(1)
            stat = popen('qstat|grep '+str(workid)).read().split()
            if len(stat) >= 4 and stat[4].find('E') != -1:
                popen('qdel '+str(workid))
    else:
        core = core if core else Getpronum()
        temp = popen('mpicc '+path.join(_path__, '..', '.Cache', jobnumber+'.c')+' -O3 '
                     '-lm -o '+path.join(_path__, '..', '.Cache', jobnumber+'.out')).read()
        popen('mpiexec -np '+str(core)+' '+path.join(_path__,
                                                     '..', '.Cache', jobnumber+'.out')+' 1>&2').read()


def Clearjob(jobnumber):
    '''Delete job files for a sepecific jobnumber.
    '''
    import os
    flist = os.listdir(os.path.join(_path__, '..', '.Cache'))
    for i in flist:
        if i.startswith(jobnumber+'.'):
            os.remove(os.path.join(_path__, '..', '.Cache', i))
    if os.listdir(os.path.join(_path__, '..', '.Cache')) == []:
        os.removedirs(os.path.join(_path__, '..', '.Cache'))
    flist = os.listdir(os.path.join(_path__, '..', '.Result'))
    for i in flist:
        if i.startswith(jobnumber+'-'):
            os.remove(os.path.join(_path__, '..', '.Result', i))
    if os.listdir(os.path.join(_path__, '..', '.Result')) == []:
        os.removedirs(os.path.join(_path__, '..', '.Result'))


_Startup()
