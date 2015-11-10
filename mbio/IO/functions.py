# -*- coding: utf-8 -*-
"""This module contains some general IO functions.
"""

__author__ = 'Wenzhi Mao'

__all__ = ["parseMRC", "parseMRCHeader", "writeMRC", "parsePDB",
           "writePDB"]


def parseMRC(filename=None, **kwargs):
    """Parse the MRC from a file."""

    from .mrc import MRC
    from os.path import exists, isfile

    if type(filename) == type(None):
        from .output import printError
        printError("The filename is wrong.")
        return None
    if exists(filename) and isfile(filename):
        return MRC(filename=filename)
    else:
        from .output import printError
        printError("The filename doesn't exists or is not a file.")


def parseMRCHeader(filename=None, **kwargs):
    """Parse the MRC header from a file."""

    from .mrc import MRCHeader

    from os.path import exists, isfile

    if type(filename) == type(None):
        from .output import printError
        printError("The filename is wrong.")
        return None
    if exists(filename) and isfile(filename):
        return MRCHeader(filename=filename)
    else:
        from .output import printError
        printError("The filename doesn't exists or is not a file.")


def writeMRC(filename, mrc, **kwargs):
    """Write the MRC file into a file."""

    from .mrc import MRC

    if not isinstance(mrc, MRC):
        from .output import printError
        printError("The mrc is not a mbio MRC instance.")
        return None
    mrc.writeData(filename=filename, **kwargs)


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
        mapindex = [i for i in range(len(names)) if not names[i].lower().endswith(
            ".pdb") and not names[i].lower().endswith('.pdb.beta')]
        pdbs = [
            i for i in range(len(names)) if names[i].lower().endswith('.pdb')]
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
            for i in range(len(pdbtemp)):
                pdbtemp[i].setChids([mapper[mapperkeylist[i]][j] if j in mapper[
                                    mapperkeylist[i]].keys() else j for j in pdbtemp[i].getChids()])
            pdb = pdbtemp[0]
            for i in range(1, len(pdbtemp)):
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
    from numpy import array
    from os.path import exists, isfile
    from ..IO.output import printInfo

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
                    'Write the beta values to {0}'.format(filename + ".beta"))
            except:
                pass
