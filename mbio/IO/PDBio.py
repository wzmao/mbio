# -*- coding: utf-8 -*-
"""This module contains some PDB IO functions.
"""

__all__ = ["parsePDB", "writePDB"]


def parsePDB(filename, **kwargs):
    """Read PDB file. If there exists a beta file read it also.

    If the file is a `.tar.gz` or `.tar.bz2` file, the file will be treated as PDB bundle file.
    Any non-PDB file will be read as the chain mapper.
        If no chain mapper then only 1 PDB file allowed in the tar file.
        If more than 1 mapper found, error raised.
        1 mapper file found, Read the mapper and read PDB files. The chain will be changed as mapper said.
    If beta file exists in the tar file, it will also be read.
    """

    from prody import parsePDB as p
    from prody import parsePDBStream as ps
    from os.path import exists, isfile, splitext
    import tarfile
    from ..IO.output import printError

    if not filename.lower().endswith('tar.gz') and not filename.lower().endswith('tar.bz2'):
        pdb = p(filename, **kwargs)
        if exists(filename + ".beta") and isfile(filename + ".beta"):
            try:
                f = open(filename + ".beta", 'r')
                beta = [float(i) for i in f]
                f.close()
                pdb.setBetas(beta)
            except:
                pass
    else:
        t = tarfile.open(
            filename, 'r:{0}'.format(filename.lower().split('.')[-1]))
        files = [i for i in t.getmembers() if i.isfile()]
        names = [i.name for i in files]
        mapindex = [i for i in xrange(len(names)) if not names[i].lower().endswith(
            ".pdb") and not names[i].lower().endswith('.pdb.beta')]
        pdbs = [
            i for i in xrange(len(names)) if names[i].lower().endswith('.pdb')]
        if len(mapindex) > 1:
            printError(
                "There are more than 1 mapper file, please check the tar.")
            return None
        elif len(mapindex) == 0:
            if len(pdbs) != 1:
                printError(
                    "No mapper and more than 1 pdb file, please check the tar.")
                return None
            else:
                title = kwargs.get('title', None)
                if title is None:
                    title, ext = splitext(names[pdbs[0]])
                    if len(title) == 7 and title.startswith('pdb'):
                        title = title[3:]
                    kwargs['title'] = title
                pdb = ps(t.extractfile(files[pdbs[0]]), **kwargs)
                if names[pdbs[0]] + '.beta' in names:
                    try:
                        f = t.extractfile(
                            files[names.index(names[pdbs[0]] + '.beta')])
                        beta = [float(i) for i in f]
                        f.close()
                        pdb.setBetas(beta)
                    except:
                        printError(
                            'The Beta file of {0} could not be parsed.'.format(names[pdbs[0]]))
        elif len(mapindex) == 1:
            mapindex = mapindex[0]
            f = t.extractfile(files[mapindex])
            nowfile = None
            mapper = {}
            for i in f:
                if i.find('New chain ID') != -1:
                    continue
                i = i.strip()
                if i:
                    if i.endswith(":"):
                        nowfile = i.replace(':', '')
                        if not nowfile in names:
                            printError(
                                "{0} is in the mapper file, but does not exists.".format(nowfile))
                            return None
                        if not nowfile in mapper.keys():
                            mapper[nowfile] = {}
                    elif len(i.split()) == 2:
                        mapper[nowfile][i.split()[0]] = i.split()[1]
            f.close()
            pdbtemp = []
            mapperkeylist = mapper.keys()
            mapperkeylist.sort()
            for mapperfile in mapperkeylist:
                f = t.extractfile(files[names.index(mapperfile)])
                pdbtemp.append(ps(f, **kwargs))
                f.close()
                if names[names.index(mapperfile)] + '.beta' in names:
                    try:
                        f = t.extractfile(
                            files[names.index(names[names.index(mapperfile)] + '.beta')])
                        beta = [float(i) for i in f]
                        f.close()
                        pdbtemp[-1].setBetas(beta)
                    except:
                        printError(
                            'The Beta file of {0} could not be parsed.'.format(names.index(mapperfile)))
            if len(set([i.numCoordsets() for i in pdbtemp])) != 1:
                printError(
                    "Each PDB must provide the same number of coordinate sets.")
                return None
            for i in xrange(len(pdbtemp)):
                pdbtemp[i].setChids([mapper[mapperkeylist[i]][j] if j in mapper[
                                    mapperkeylist[i]].keys() else j for j in pdbtemp[i].getChids()])
            pdb = pdbtemp[0]
            for i in xrange(1, len(pdbtemp)):
                pdb = pdb + pdbtemp[i]
            if filename.lower().endswith('.tar.gz'):
                filename = filename[:-7]
            elif filename.lower().endswith('.tar.bz2'):
                filename = filename[:-8]
            pdb.setTitle(filename)
        t.close()
    return pdb


def writePDB(filename, pdb, forcebeta=None, **kwargs):
    """Write PDB file. Write the beta to beta also."""

    from prody import writePDB as w
    from prody import writePDBStream as ws
    import tarfile
    from numpy import array, dtype
    from os.path import exists, isfile, split
    from ..IO.output import printInfo, printError
    import StringIO

    if set([len(i) for i in set(pdb.getChids())]) != set([1]) or pdb.numAtoms() > 99999:
        gotar = True
    else:
        if filename.lower().endswith('.pdb') or filename.lower().endswith('.pdb.gz'):
            gotar = False
        else:
            gotar = True

    if gotar:
        if filename.lower().endswith('.tar.gz') or filename.lower().endswith('.tar.bz2'):
            pass
        elif filename.lower().endswith('.pdb'):
            filename = filename + '.tar.gz'
        elif filename.lower().endswith('.gz') and not filename.lower().endswith('.tar.gz'):
            filename = '.'.join(
                filename.split('.')[:-1] + ['tar'] + filename.split('.')[-1:])
        elif filename.lower().endswith('.bz2') and not filename.lower().endswith('.tar.bz2'):
            filename = '.'.join(
                filename.split('.')[:-1] + ['tar'] + filename.split('.')[-1:])
        elif filename.lower().endswith('.tar'):
            filename = filename + '.gz'
        else:
            filename = filename + '.tar.gz'

    if filename.lower().endswith('.tar.gz') or filename.lower().endswith('.tar.bz2'):
        filetype = filename.lower().split('.')[-1]
        rootname = split(filename)[1]
        rootname = rootname[:-7] if filetype == 'gz' else rootname[:-8]
        chains = list(set(pdb.getChids()))
        chains.sort()
        chainnums = [pdb.select('chain ' + i).numAtoms() for i in chains]
        sumatom = 0
        segs = [[]]
        for i in chains:
            segs[-1].append(i)
            sumatom += chainnums[chains.index(i)]
            onetoobig = False
            if (sumatom > 99999) or len(segs[-1]) > 26:
                if len(segs[-1]) != 1:
                    segs.append([segs[-1][-1]])
                    segs[-2] = segs[-2][:-1]
                    sumatom = chainnums[chains.index(i)]
                    if (chainnums[chains.index(i)] > 99999):
                        onetoobig = True
                else:
                    onetoobig = True
            if onetoobig or len(segs[-1]) == 26:
                sumatom = 0
                segs.append([])
        while segs[-1] == []:
            segs = segs[:-1]
        t = tarfile.open(filename, "w:" + filetype)
        nownum = 1
        info = t.tarinfo()
        mapper = {}

        for i in xrange(len(segs)):
            p1 = pdb.select('chain ' + ' '.join(segs[i]))
            if p1.numAtoms() < 99999:
                nowfilename = rootname + '-bundle-' + str(nownum) + '.pdb'
                mapper[nowfilename] = {}
                for j in xrange(len(segs[i])):
                    mapper[nowfilename][segs[i][j]] = chr(ord('A') + j)
                p1.setChids([mapper[nowfilename][j] for j in p1.getChids()])

                data = StringIO.StringIO()
                ws(data, p1)
                data.seek(0)

                info.name = nowfilename
                info.size = data.len
                info.mtime = tarfile.time.time()

                t.addfile(info, data)
                nownum += 1

                if forcebeta is None:
                    if (array([float('%6.2f' % i) for i in p1.getBetas()]) == p1.getBetas()).all():
                        forcebeta1 = False
                    else:
                        forcebeta1 = True
                else:
                    forcebeta1 = forcebeta
                if forcebeta1:
                    try:
                        data = StringIO.StringIO()
                        data.write('\n'.join([str(i) for i in p1.getBetas()]))
                        data.seek(0)
                        info.name = nowfilename + '.beta'
                        info.size = data.len
                        info.mtime = tarfile.time.time()
                        t.addfile(info, data)
                        printInfo(
                            'Write the Beta values to {0}'.format(nowfilename + ".beta"))
                    except:
                        printError(
                            "Writting the Beta for {0} failed.".format(nowfilename))
                        pass
            else:
                ss = 0
                while p1[ss:ss + 99999]:
                    nowfilename = rootname + '-bundle-' + str(nownum) + '.pdb'
                    mapper[nowfilename] = {}
                    mapper[nowfilename][segs[i][0]] = 'A'
                    p2 = p1[ss:ss + 99999]
                    p2.setChids([mapper[nowfilename][j]
                                 for j in p2.getChids()])

                    data = StringIO.StringIO()
                    ws(data, p2)
                    data.seek(0)

                    info.name = nowfilename
                    info.size = data.len
                    info.mtime = tarfile.time.time()

                    t.addfile(info, data)
                    nownum += 1
                    ss += 99999
                    if forcebeta is None:
                        if (array([float('%6.2f' % i) for i in p2.getBetas()]) == p2.getBetas()).all():
                            forcebeta1 = False
                        else:
                            forcebeta1 = True
                    else:
                        forcebeta1 = forcebeta
                    if forcebeta1:
                        try:
                            data = StringIO.StringIO()
                            data.write(
                                '\n'.join([str(i) for i in p2.getBetas()]))
                            data.seek(0)
                            info.name = nowfilename + '.beta'
                            info.size = data.len
                            info.mtime = tarfile.time.time()
                            t.addfile(info, data)
                            printInfo(
                                'Write the Beta values to {0}'.format(nowfilename + ".beta"))
                        except:
                            printError(
                                "Writting the Beta for {0} failed.".format(nowfilename))
                            pass

        data = StringIO.StringIO()
        data.write("    New chain ID            Original chain ID\n")
        for i in sorted(mapper.keys()):
            data.write('\n' + i + ':\n')
            for j in sorted(mapper[i].keys()):
                data.write(
                    "           {0}           {1}\n".format(j, mapper[i][j]))
        data.seek(0)
        info.name = rootname + "-chain-id-mapping.txt"
        info.size = data.len
        info.mtime = tarfile.time.time()
        t.addfile(info, data)
        t.close()
        return filename
    else:
        w(filename, pdb, **kwargs)
        if forcebeta is None:
            if (array([float('%6.2f' % i) for i in pdb.getBetas()]) == pdb.getBetas()).all():
                forcebeta = False
            else:
                forcebeta = True
        if forcebeta:
            if (exists(filename + ".beta") and isfile(filename + ".beta")) or not exists(filename + '.beta'):
                try:
                    f = open(filename + ".beta", 'w')
                    f.write('\n'.join([str(i) for i in pdb.getBetas()]))
                    f.close()
                    printInfo(
                        'Write the Beta values to {0}'.format(filename + ".beta"))
                except:
                    printError(
                        "Writting the Beta for {0} failed.".format(filename))
                    pass
        return filename
