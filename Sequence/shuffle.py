'''OMES is a simple protein correlation algrithm.
We use OMES to calculate the correlation for 20 kinds of amino acid and gaps.(21)
It also provide a shuffling function to calculate the shuffled P-value.
The shuffled P-value calculation is performced by C and MPI.
'''
__author__ = 'Wenzhi Mao'
__all__ = ['ShuffleOMES', 'ShuffleMI', 'ShuffleMIp']


def _Startup():
    from mbio import _ABSpath
    global _path__
    _path__ = _ABSpath()


def ShuffleOMES(sequences, times=10000, cutoff=0.2, core=1, output=1, cluster=0, save=False):
    '''It is a function to calculate the p value for shuffled OMES.
    Given the sequences in a list with no format.
    times is the shuffle times.
    cutoff is the lower cutoff. (Haven't finished yet. Calc all now.)
    core is the process number to run.
    output is a mark for output.'''
    import os
    from os import path
    import ctypes as ct
    import struct
    if not os.path.exists(_path__+'/../.Cache'):
        os.mkdir(_path__+'/../.Cache')
    plist = os.listdir(_path__+'/../.Cache')
    i = 1
    while os.path.exists(_path__+'/../.Cache/'+str(i)+'.fasta'):
        i += 1
    f = open(_path__+'/../.Cache/'+str(i)+'.fasta', 'w')
    fastatempname = str(i)+'.fasta'
    f.write('\n'.join(sequences))
    f.close()
    if not os.path.exists(_path__+'/../.Result'):
        os.mkdir(_path__+'/../.Result')
    f = open(_path__+'/omes_mpi.c')
    script = f.read()
    f.close()
    if output:
        output = 1
    else:
        output = 0
    f = open(_path__+'/../.Cache/'+str(i)+'.c', 'w')
    f.write(script.replace
            ('#define seqnum 200', '#define seqnum '+str(len(sequences)))
            .replace
            ('#define lennum 700', '#define lennum '+str(len(sequences[0])+1))
            .replace
            ('file.fasta', _path__+'/../.Cache/'+str(i)+'.fasta')
            .replace
            ('OUTPUT 1', 'OUTPUT '+str(output))
            .replace
            ('#define times 10000', '#define times '+str(times))
            .replace
            ('#define cutoff 0.8', '#define cutoff '+str(1-cutoff))
            .replace
            ("OMESsave.save", _path__+'/../.Result/'+str(i)+'-omes.save')
            .replace
            ("Psave.save", _path__+'/../.Result/'+str(i)+'-p.save'))
    f.close()
    if not cluster:
        temp = os.popen('mpicc '+_path__+'/../.Cache/'+str(i)+'.c -O3 '
                        '-lm -o '+_path__+'/../.Cache/'+str(i)+'.out').read()
        os.popen('mpiexec -np '+str(core)+' '+_path__+'/../.Cache/'+str(i) +
                 '.out 1>&2').read()
    else:
        if os.popen('hostname\n').read() == 'n000\n':
            f = open(_path__+'/../Scripts/qsub.job', 'r')
            qsub = f.read()
            f.close()
            find = False
            while not find:
                for core in [64, 32, 24, 12, 8, 4]:
                    if int(os.popen('qstat -f | grep BP|grep parallel.q |grep "/{0} "|grep "/0/{0} "|grep -v "au"|grep -v "o"|grep "" -c'.format(str(core))).read()) != 0:
                        find = True
                        break
            clusternames = ['bahar64', 'core32', 'core24', 'bahar12',
                            'core8', 'core4']
            clustername = clusternames[[64, 32, 24, 12, 8, 4].index(core)]
            qsub = [qsub.replace('core_number', str(core)).replace
                    ('clustername', clustername).replace
                    ('output_name', str(i)+'.output').replace
                    ('error_name', str(i)+'.err').replace
                    ('c_file_name', str(i)+'.c').replace
                    ('out_file_name', str(i)+'.out').replace
                    ('job_name', 'job_'+str(i))][0]
            f = open(_path__+'/../.Cache/'+str(i)+'.job', 'w')
            f.write(qsub)
            f.close()
            workid = os.popen(
                'cd '+_path__+'/../.Cache/;qsub '+str(i)+'.job').read()
            workid = workid.split()[2]
            import time
            while os.popen('qstat -j '+str(workid)+' 2>&1').read().startswith('===='):
                time.sleep(1)
                stat = os.popen('qstat|grep '+str(workid)).read().split()
                if len(stat) >= 4 and stat[4].find('E') != -1:
                    os.popen('qdel '+str(workid))
        else:
            print '* The system cannot be recognized.'
            return None
    if path.exists(_path__+'/../.Result/'+str(i)+'-omes.save'):
        f = open(_path__+'/../.Result/'+str(i)+'-omes.save', 'r')
        tempom = f.read()
        f.close()
        if len(tempom) != 8*len(sequences[0])**2:
            print '* OMES file size wrong.'
            return None
        om = []
        k = 0
        for kk in range(len(sequences[0])):
            om.append([])
            for j in range(len(sequences[0])):
                om[-1].append(struct.unpack('d', tempom[k:k+8])[0])
                k += 8
        if save:
            import shutil
            shutil.copy(_path__+'/../.Result/'+str(i)+'-omes.save',  './')
            os.rename('./'+str(i)+'-omes.save', './omes.save')
        os.remove(_path__+'/../.Result/'+str(i)+'-omes.save')
    else:
        return None
    if path.exists(_path__+'/../.Result/'+str(i)+'-p.save'):
        f = open(_path__+'/../.Result/'+str(i)+'-p.save', 'r')
        tempp = f.read()
        f.close()
        if len(tempp) != 4 * len(sequences[0])**2:
            print '* P file size wrong.'
            return None
        p = []
        k = 0
        for kk in range(len(sequences[0])):
            p.append([])
            for j in range(len(sequences[0])):
                p[-1].append(struct.unpack('i', tempp[k:k+4])[0])
                k += 4
        if save:
            import shutil
            shutil.copy(_path__+'/../.Result/'+str(i)+'-p.save',  './')
            os.rename('./'+str(i)+'-p.save', './p.save')
        os.remove(_path__+'/../.Result/'+str(i)+'-p.save')
    else:
        return None
    if path.exists(_path__+'/../.Cache/'+str(i)+'.fasta'):
        os.remove(_path__+'/../.Cache/'+str(i)+'.fasta')
    if path.exists(_path__+'/../.Cache/'+str(i)+'.c'):
        os.remove(_path__+'/../.Cache/'+str(i)+'.c')
    if path.exists(_path__+'/../.Cache/'+str(i)+'.out'):
        os.remove(_path__+'/../.Cache/'+str(i)+'.out')
    if len(os.listdir(_path__+'/../.Cache')) == 0:
        os.removedirs(_path__+'/../.Cache/')
        os.removedirs(_path__+'/../.Result/')
    return p


def ShuffleMI(sequences, times=10000, cutoff=0.2, core=1, output=1, cluster=False, save=False):
    '''It is a function to calculate the p value for shuffled MI.
    Given the sequences in a list with no format.
    times is the shuffle times.
    cutoff is the lower cutoff. (Haven't finished yet. Calc all now.)
    core is the process number to run.
    output is a mark for output.'''
    import os
    from os import path
    import ctypes as ct
    import struct
    if not os.path.exists(_path__+'/../.Cache'):
        os.mkdir(_path__+'/../.Cache')
    plist = os.listdir(_path__+'/../.Cache')
    i = 1
    while os.path.exists(_path__+'/../.Cache/'+str(i)+'.fasta'):
        i += 1
    f = open(_path__+'/../.Cache/'+str(i)+'.fasta', 'w')
    fastatempname = str(i)+'.fasta'
    f.write('\n'.join(sequences))
    f.close()
    if not os.path.exists(_path__+'/../.Result'):
        os.mkdir(_path__+'/../.Result')
    f = open(_path__+'/mi_mpi.c')
    script = f.read()
    f.close()
    if output:
        output = 1
    else:
        output = 0
    f = open(_path__+'/../.Cache/'+str(i)+'.c', 'w')
    f.write(script.replace
            ('#define seqnum 200', '#define seqnum '+str(len(sequences)))
            .replace
            ('#define lennum 700', '#define lennum '+str(len(sequences[0])+1))
            .replace
            ('file.fasta', _path__+'/../.Cache/'+str(i)+'.fasta')
            .replace
            ('OUTPUT 1', 'OUTPUT '+str(output))
            .replace
            ('#define times 10000', '#define times '+str(times))
            .replace
            ('#define cutoff 0.8', '#define cutoff '+str(1-cutoff))
            .replace
            ("MIsave.save", _path__+'/../.Result/'+str(i)+'-mi.save')
            .replace
            ("Psave.save", _path__+'/../.Result/'+str(i)+'-p.save'))
    f.close()
    if not cluster:
        temp = os.popen('mpicc '+_path__+'/../.Cache/'+str(i)+'.c -O3 '
                        '-lm -o '+_path__+'/../.Cache/'+str(i)+'.out').read()
        os.popen('mpiexec -np '+str(core)+' '+_path__+'/../.Cache/'+str(i) +
                 '.out 1>&2').read()
    else:
        if os.popen('hostname\n').read() == 'n000\n':
            f = open(_path__+'/../Scripts/qsub.job', 'r')
            qsub = f.read()
            f.close()
            find = False
            while not find:
                for core in [64, 32, 24, 12, 8, 4]:
                    if int(os.popen('qstat -f | grep BP|grep parallel.q |grep "/{0} "|grep "/0/{0} "|grep -v "au"|grep -v "o"|grep "" -c'.format(str(core))).read()) != 0:
                        find = True
                        break
            clusternames = ['bahar64', 'core32', 'core24', 'bahar12',
                            'core8', 'core4']
            clustername = clusternames[[64, 32, 24, 12, 8, 4].index(core)]
            qsub = [qsub.replace('core_number', str(core)).replace
                    ('clustername', clustername).replace
                    ('output_name', str(i)+'.output').replace
                    ('error_name', str(i)+'.err').replace
                    ('c_file_name', str(i)+'.c').replace
                    ('out_file_name', str(i)+'.out').replace
                    ('job_name', 'job_'+str(i))][0]
            f = open(_path__+'/../.Cache/'+str(i)+'.job', 'w')
            f.write(qsub)
            f.close()
            workid = os.popen(
                'cd '+_path__+'/../.Cache/;qsub '+str(i)+'.job').read()
            workid = workid.split()[2]
            import time
            while os.popen('qstat -j '+str(workid)+' 2>&1').read().startswith('===='):
                time.sleep(1)
                stat = os.popen('qstat|grep '+str(workid)).read().split()
                if len(stat) >= 4 and stat[4].find('E') != -1:
                    os.popen('qdel '+str(workid))
        else:
            print '* The system cannot be recognized.'
            return None
    if path.exists(_path__+'/../.Result/'+str(i)+'-mi.save'):
        f = open(_path__+'/../.Result/'+str(i)+'-mi.save', 'r')
        tempmi = f.read()
        f.close()
        if len(tempmi) != 8*len(sequences[0])**2:
            print '* Mi file size wrong.'
            return None
        mi = []
        k = 0
        tempmi = struct.unpack(str(int(len(tempmi)/8))+'d', tempmi)
        l = len(sequences[0])
        for ii in range(l):
            mi.append(tempmi[ii*l:ii*l+l])
        if save:
            import shutil
            shutil.copy(_path__+'/../.Result/'+str(i)+'-mi.save',  './')
            os.rename('./'+str(i)+'-mi.save', './mi.save')
        os.remove(_path__+'/../.Result/'+str(i)+'-mi.save')
    else:
        return None
    if path.exists(_path__+'/../.Result/'+str(i)+'-p.save'):
        f = open(_path__+'/../.Result/'+str(i)+'-p.save', 'r')
        tempp = f.read()
        f.close()
        if len(tempp) != 4 * len(sequences[0])**2:
            print '* P file size wrong.'
            return None
        p = []
        k = 0
        tempp = struct.unpack(str(int(len(tempp)/4))+'i', tempp)
        l = len(sequences[0])
        for ii in range(l):
            p.append(tempp[ii*l:ii*l+l])
        if save:
            import shutil
            shutil.copy(_path__+'/../.Result/'+str(i)+'-p.save',  './')
            os.rename('./'+str(i)+'-p.save', './p.save')
        os.remove(_path__+'/../.Result/'+str(i)+'-p.save')
    else:
        return None
    if path.exists(_path__+'/../.Cache/'+str(i)+'.fasta'):
        os.remove(_path__+'/../.Cache/'+str(i)+'.fasta')
    if path.exists(_path__+'/../.Cache/'+str(i)+'.c'):
        os.remove(_path__+'/../.Cache/'+str(i)+'.c')
    if path.exists(_path__+'/../.Cache/'+str(i)+'.out'):
        os.remove(_path__+'/../.Cache/'+str(i)+'.out')
    if path.exists(_path__+'/../.Cache/'+str(i)+'.output'):
        os.remove(_path__+'/../.Cache/'+str(i)+'.output')
    if path.exists(_path__+'/../.Cache/'+str(i)+'.job'):
        os.remove(_path__+'/../.Cache/'+str(i)+'.job')
    if path.exists(_path__+'/../.Cache/'+str(i)+'.err'):
        os.remove(_path__+'/../.Cache/'+str(i)+'.err')
    if len(os.listdir(_path__+'/../.Cache')) == 0:
        os.removedirs(_path__+'/../.Cache/')
        os.removedirs(_path__+'/../.Result/')
    return p


def ShuffleMIp(sequences, times=10000, cutoff=0.2, core=1, output=1, cluster=False, save=False):
    '''It is a function to calculate the p value for shuffled MIp.
    Given the sequences in a list with no format.
    times is the shuffle times.
    cutoff is the lower cutoff. (Haven't finished yet. Calc all now.)
    core is the process number to run.
    output is a mark for output.'''
    import os
    from os import path
    import ctypes as ct
    import struct
    if not os.path.exists(_path__+'/../.Cache'):
        os.mkdir(_path__+'/../.Cache')
    plist = os.listdir(_path__+'/../.Cache')
    i = 1
    while os.path.exists(_path__+'/../.Cache/'+str(i)+'.fasta'):
        i += 1
    f = open(_path__+'/../.Cache/'+str(i)+'.fasta', 'w')
    fastatempname = str(i)+'.fasta'
    f.write('\n'.join(sequences))
    f.close()
    if not os.path.exists(_path__+'/../.Result'):
        os.mkdir(_path__+'/../.Result')
    f = open(_path__+'/mip_mpi.c')
    script = f.read()
    f.close()
    if output:
        output = 1
    else:
        output = 0
    f = open(_path__+'/../.Cache/'+str(i)+'.c', 'w')
    f.write(script.replace
            ('#define seqnum 200', '#define seqnum '+str(len(sequences)))
            .replace
            ('#define lennum 700', '#define lennum '+str(len(sequences[0])+1))
            .replace
            ('file.fasta', _path__+'/../.Cache/'+str(i)+'.fasta')
            .replace
            ('OUTPUT 1', 'OUTPUT '+str(output))
            .replace
            ('#define times 10000', '#define times '+str(times))
            .replace
            ('#define cutoff 0.8', '#define cutoff '+str(1-cutoff))
            .replace
            ("MIpsave.save", _path__+'/../.Result/'+str(i)+'-mip.save')
            .replace
            ("Psave.save", _path__+'/../.Result/'+str(i)+'-p.save'))
    f.close()
    if not cluster:
        temp = os.popen('mpicc '+_path__+'/../.Cache/'+str(i)+'.c -O3 '
                        '-lm -o '+_path__+'/../.Cache/'+str(i)+'.out').read()
        os.popen('mpiexec -np '+str(core)+' '+_path__+'/../.Cache/'+str(i) +
                 '.out 1>&2').read()
    else:
        if os.popen('hostname\n').read() == 'n000\n':
            f = open(_path__+'/../Scripts/qsub.job', 'r')
            qsub = f.read()
            f.close()
            find = False
            while not find:
                for core in [64, 32, 24, 12, 8, 4]:
                    if int(os.popen('qstat -f | grep BP|grep parallel.q |grep "/{0} "|grep "/0/{0} "|grep -v "au"|grep -v "o"|grep "" -c'.format(str(core))).read()) != 0:
                        find = True
                        break
            clusternames = ['bahar64', 'core32', 'core24', 'bahar12',
                            'core8', 'core4']
            clustername = clusternames[[64, 32, 24, 12, 8, 4].index(core)]
            qsub = [qsub.replace('core_number', str(core)).replace
                    ('clustername', clustername).replace
                    ('output_name', str(i)+'.output').replace
                    ('error_name', str(i)+'.err').replace
                    ('c_file_name', str(i)+'.c').replace
                    ('out_file_name', str(i)+'.out').replace
                    ('job_name', 'job_'+str(i))][0]
            f = open(_path__+'/../.Cache/'+str(i)+'.job', 'w')
            f.write(qsub)
            f.close()
            workid = os.popen(
                'cd '+_path__+'/../.Cache/;qsub '+str(i)+'.job').read()
            workid = workid.split()[2]
            import time
            while os.popen('qstat -j '+str(workid)+' 2>&1').read().startswith('===='):
                time.sleep(1)
                stat = os.popen('qstat|grep '+str(workid)).read().split()
                if len(stat) >= 4 and stat[4].find('E') != -1:
                    os.popen('qdel '+str(workid))
        else:
            print '* The system cannot be recognized.'
            return None
    if path.exists(_path__+'/../.Result/'+str(i)+'-mip.save'):
        f = open(_path__+'/../.Result/'+str(i)+'-mip.save', 'r')
        tempmi = f.read()
        f.close()
        if len(tempmi) != 8*len(sequences[0])**2:
            print '* Mip file size wrong.'
            return None
        mi = []
        k = 0
        tempmi = struct.unpack(str(int(len(tempmi)/8))+'d', tempmi)
        l = len(sequences[0])
        for ii in range(l):
            mi.append(tempmi[ii*l:ii*l+l])
        if save:
            import shutil
            shutil.copy(_path__+'/../.Result/'+str(i)+'-mip.save',  './')
            os.rename('./'+str(i)+'-mip.save', './mip.save')
        os.remove(_path__+'/../.Result/'+str(i)+'-mip.save')
    else:
        return None
    if path.exists(_path__+'/../.Result/'+str(i)+'-p.save'):
        f = open(_path__+'/../.Result/'+str(i)+'-p.save', 'r')
        tempp = f.read()
        f.close()
        if len(tempp) != 4 * len(sequences[0])**2:
            print '* P file size wrong.'
            return None
        p = []
        k = 0
        tempp = struct.unpack(str(int(len(tempp)/4))+'i', tempp)
        l = len(sequences[0])
        for ii in range(l):
            p.append(tempp[ii*l:ii*l+l])
        if save:
            import shutil
            shutil.copy(_path__+'/../.Result/'+str(i)+'-p.save',  './')
            os.rename('./'+str(i)+'-p.save', './p.save')
        os.remove(_path__+'/../.Result/'+str(i)+'-p.save')
    else:
        return None
    if path.exists(_path__+'/../.Cache/'+str(i)+'.fasta'):
        os.remove(_path__+'/../.Cache/'+str(i)+'.fasta')
    if path.exists(_path__+'/../.Cache/'+str(i)+'.c'):
        os.remove(_path__+'/../.Cache/'+str(i)+'.c')
    if path.exists(_path__+'/../.Cache/'+str(i)+'.out'):
        os.remove(_path__+'/../.Cache/'+str(i)+'.out')
    if path.exists(_path__+'/../.Cache/'+str(i)+'.output'):
        os.remove(_path__+'/../.Cache/'+str(i)+'.output')
    if path.exists(_path__+'/../.Cache/'+str(i)+'.job'):
        os.remove(_path__+'/../.Cache/'+str(i)+'.job')
    if path.exists(_path__+'/../.Cache/'+str(i)+'.err'):
        os.remove(_path__+'/../.Cache/'+str(i)+'.err')
    if len(os.listdir(_path__+'/../.Cache')) == 0:
        os.removedirs(_path__+'/../.Cache/')
        os.removedirs(_path__+'/../.Result/')
    return p


_Startup()