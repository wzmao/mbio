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


def ShuffleOMES(sequences, times=10000, cutoff=0.2, core=0, output=1, cluster=0, save=False):
    '''It is a function to calculate the p value for shuffled OMES.
    Given the sequences in a list with no format.
    times is the shuffle times.
    cutoff is the lower cutoff. (Haven't finished yet. Calc all now.)
    core is the process number to run.
    output is a mark for output.'''
    import os
    from os import path
    from mbio.Application import job_organization as jo
    scriptfile=path.join(_path__,'..','Scripts','omes_mpi.c')
    jobnumber=jo.AskJobNumber()
    f = open(path.join(_path__,'..','.Cache',jobnumber+'.fasta'), 'w')
    f.write('\n'.join(sequences))
    f.close()
    jo.MkdirResult()
    f = open(scriptfile,'r')
    script = f.read()
    f.close()
    output='1' if output else '0'
    jo.Writejob(jobnumber,script.replace
            ('#define seqnum', '#define seqnum '+str(len(sequences)))
            .replace
            ('#define lennum', '#define lennum '+str(len(sequences[0])+1))
            .replace
            ('file.fasta', path.join(_path__,'..','.Cache',jobnumber+'.fasta'))
            .replace
            ('#define OUTPUT', '#define OUTPUT '+str(output))
            .replace
            ('#define times', '#define times '+str(times))
            .replace
            ('#define cutoff', '#define cutoff '+str(1-cutoff))
            .replace
            ("OMESsave.save", path.join(_path__,'..','.Result',jobnumber+'-omes.save'))
            .replace
            ("Psave.save", path.join(_path__,'..','.Result',jobnumber+'-p.save')))
    jo.SubShufflejob(jobnumber,cluster,core)
    if path.exists(path.join(_path__,'..','.Result',jobnumber+'-omes.save')):
        from mbio.IO import matrix
        m = matrix.ReadMatrix(path.join(_path__,'..','.Result',jobnumber+'-omes.save'), 'd',len(sequences[0]))
        if save:
            import shutil
            shutil.copy(path.join(_path__,'..','.Result',jobnumber+'-omes.save'),  './')
            os.rename(path.join('.',jobnumber+'-omes.save'), './omes.save')
    else:
        return None
    if path.exists(path.join(_path__,'..','.Result',jobnumber+'-p.save')):
        from mbio.IO import matrix
        p = matrix.ReadMatrix(path.join(_path__,'..','.Result',jobnumber+'-p.save'), 'i',len(sequences[0]))
        if save:
            import shutil
            shutil.copy(path.join(_path__,'..','.Result',jobnumber+'-p.save'),  './')
            os.rename(path.join('.',jobnumber+'-p.save'), './p.save')
    else:
        return None
    jo.Clearjob(jobnumber)
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
    from mbio.Application import job_organization as jo
    scriptfile=path.join(_path__,'..','Scripts','mi_mpi.c')
    jobnumber=jo.AskJobNumber()
    f = open(path.join(_path__,'..','.Cache',jobnumber+'.fasta'), 'w')
    f.write('\n'.join(sequences))
    f.close()
    jo.MkdirResult()
    f = open(scriptfile,'r')
    script = f.read()
    f.close()
    output='1' if output else '0'
    jo.Writejob(jobnumber,script.replace
            ('#define seqnum', '#define seqnum '+str(len(sequences)))
            .replace
            ('#define lennum', '#define lennum '+str(len(sequences[0])+1))
            .replace
            ('file.fasta', path.join(_path__,'..','.Cache',jobnumber+'.fasta'))
            .replace
            ('#define OUTPUT', '#define OUTPUT '+str(output))
            .replace
            ('#define times', '#define times '+str(times))
            .replace
            ('#define cutoff', '#define cutoff '+str(1-cutoff))
            .replace
            ("MIsave.save", path.join(_path__,'..','.Result',jobnumber+'-mi.save'))
            .replace
            ("Psave.save", path.join(_path__,'..','.Result',jobnumber+'-p.save')))
    jo.SubShufflejob(jobnumber,cluster,core)
    if path.exists(path.join(_path__,'..','.Result',jobnumber+'-mi.save')):
        from mbio.IO import matrix
        m = matrix.ReadMatrix(path.join(_path__,'..','.Result',jobnumber+'-mi.save'), 'd',len(sequences[0]))
        if save:
            import shutil
            shutil.copy(path.join(_path__,'..','.Result',jobnumber+'-mi.save'),  './')
            os.rename(path.join('.',jobnumber+'-mi.save'), './mi.save')
    else:
        return None
    if path.exists(path.join(_path__,'..','.Result',jobnumber+'-p.save')):
        from mbio.IO import matrix
        p = matrix.ReadMatrix(path.join(_path__,'..','.Result',jobnumber+'-p.save'), 'i',len(sequences[0]))
        if save:
            import shutil
            shutil.copy(path.join(_path__,'..','.Result',jobnumber+'-p.save'),  './')
            os.rename(path.join('.',jobnumber+'-p.save'), './p.save')
    else:
        return None
    jo.Clearjob(jobnumber)
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
    from mbio.Application import job_organization as jo
    scriptfile=path.join(_path__,'..','Scripts','mip_mpi.c')
    jobnumber=jo.AskJobNumber()
    f = open(path.join(_path__,'..','.Cache',jobnumber+'.fasta'), 'w')
    f.write('\n'.join(sequences))
    f.close()
    jo.MkdirResult()
    f = open(scriptfile,'r')
    script = f.read()
    f.close()
    output='1' if output else '0'
    jo.Writejob(jobnumber,script.replace
            ('#define seqnum', '#define seqnum '+str(len(sequences)))
            .replace
            ('#define lennum', '#define lennum '+str(len(sequences[0])+1))
            .replace
            ('file.fasta', path.join(_path__,'..','.Cache',jobnumber+'.fasta'))
            .replace
            ('#define OUTPUT', '#define OUTPUT '+str(output))
            .replace
            ('#define times', '#define times '+str(times))
            .replace
            ('#define cutoff', '#define cutoff '+str(1-cutoff))
            .replace
            ("MIpsave.save", path.join(_path__,'..','.Result',jobnumber+'-mip.save'))
            .replace
            ("Psave.save", path.join(_path__,'..','.Result',jobnumber+'-p.save')))
    jo.SubShufflejob(jobnumber,cluster,core)
    if path.exists(path.join(_path__,'..','.Result',jobnumber+'-mip.save')):
        from mbio.IO import matrix
        m = matrix.ReadMatrix(path.join(_path__,'..','.Result',jobnumber+'-mip.save'), 'd',len(sequences[0]))
        if save:
            import shutil
            shutil.copy(path.join(_path__,'..','.Result',jobnumber+'-mip.save'),  './')
            os.rename(path.join('.',jobnumber+'-mip.save'), './mip.save')
    else:
        return None
    if path.exists(path.join(_path__,'..','.Result',jobnumber+'-p.save')):
        from mbio.IO import matrix
        p = matrix.ReadMatrix(path.join(_path__,'..','.Result',jobnumber+'-p.save'), 'i',len(sequences[0]))
        if save:
            import shutil
            shutil.copy(path.join(_path__,'..','.Result',jobnumber+'-p.save'),  './')
            os.rename(path.join('.',jobnumber+'-p.save'), './p.save')
    else:
        return None
    jo.Clearjob(jobnumber)
    return p


_Startup()
