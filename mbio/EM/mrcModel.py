# -*- coding: utf-8 -*-
"""This module contains some functions for EM Model.
"""

__author__ = 'Wenzhi Mao'
__all__ = ['pdb2mrc']


def GetAtomID(names):
    """Transform the atom name into id.
    H -> 0, HE -> 1 and so on"""

    from numpy import array, int32
    trans = {"H": 0, "HE": 1, "LI": 2, "BE": 3, "B": 4, "C": 5, "N": 6,
             "O": 7, "F": 8, "Ne": 9, "Na": 10, "MG": 11, "AL": 12,
             "SI": 13, "P": 14, "S": 15, "CL": 16, "AR": 17, "K": 18,
             "CA": 19, "FE": 25, }
    atomid = map(lambda x: trans.get(x[:2].strip(), 0), names)
    atomid = array(atomid, dtype=int32)
    return atomid


def pdb2mrc(pdb, grid_size, pixel_size, resolution, fsccutoff=.5, *args, **kwargs):
    """Transform PDB to MRC.
    You could pass a prody `AtomGroup` class as pdb or a path or a pdb code.

    grid_size is the number of grid point along x,y,z axis.
    pixel_size is the length between grids, the unit is Angstrom.
    resolution is the resolution of the map, the unit is Angstrom.

    You could also pass the x,y,z shift to the map, which will be added to
    the coordinate of PDB coordinates.

    The original C code is written by Xueming @ ucsf Cheng lab.
    And then optimized by Wenzhi Mao @ Tsinghua University Haipeng Gong Lab."""

    from prody import AtomGroup
    from numpy import array, argsort, zeros, float32
    from ..IO.output import printInfo, printError
    from .mrc import MRC
    from .Cmrcmodel_p import Cpdb2mrc

    if not isinstance(pdb, AtomGroup):
        from prody import parsePDB
        try:
            pdb = parsePDB(pdb, needfullname=True)
        except:
            printError("The pdb could not be recongnized. Stopped.")
            return None

    grid_size = int(grid_size)
    pixel_size = float(pixel_size)
    resolution = float(resolution)
    if len(args) >= 3:
        shift = [float(i) for i in list(args)[:3]]
        printInfo("X, Y, Z shift parsed.")
    else:
        shift = [0., 0., 0.]

    occ = pdb.getOccupancies().astype(float32)
    bf = pdb.getBetas().astype(float32)
    cor = (pdb.getCoords() + array(shift)).astype(float32)
    atomid = GetAtomID(pdb.getNames())

    if not len(occ) == len(bf) == cor.shape[0] == len(atomid):
        printError("Wrong length for all parameter. Stopped.")
        return None

    ind = argsort(atomid)
    occ = occ[ind]
    bf = bf[ind]
    cor = cor[ind]
    atomid = atomid[ind]

    mrc = MRC()
    mrc.nx = mrc.ny = mrc.nz = mrc.mx = mrc.my = mrc.mz = grid_size
    mrc.nxstart = mrc.nystart = mrc.nzstart = 0
    mrc.cella = [grid_size * pixel_size] * 3
    mrc.cellb = [90.] * 3
    mrc.mode = 2
    mrc.mapc, mrc.mapr, mrc.maps = 1, 2, 3
    mrc.dmin = mrc.dmax = mrc.dmean = 0
    mrc.ispg = 0
    mrc.nsymbt = 0
    mrc.origin = [-i for i in shift]
    mrc.map = 'MAP/0'
    mrc.machst = 0
    mrc.rms = 0.
    mrc.nlabels = 0
    mrc.label = [""] * 10
    mrc.symdata = None
    mrc.extra = ''
    mrc.xstart, mrc.ystart, mrc.zstart = array(
        mrc.origin)[argsort([mrc.mapc, mrc.mapr, mrc.maps])]
    mrc.data = zeros((grid_size, grid_size, grid_size), dtype=float32)

    if pixel_size * grid_size / resolution > grid_size // 2:
        printError("Warning: resolution is larger than Nyquist.")
        resolution = pixel_size * grid_size / (grid_size // 2)
        printError("The resolution fall back to {0}".format(resolution))

    result = Cpdb2mrc(
        atomid=atomid, occ=occ, bf=bf, cor=cor, map=mrc.data, nsam=grid_size,
        psize=pixel_size, res=resolution, fsccutoff=fsccutoff)
    if isinstance(result, (tuple, list)):
        if result[0] is None:
            printError(result[1])
            return None

    return mrc
