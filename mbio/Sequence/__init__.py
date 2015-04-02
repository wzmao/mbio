"""This module contains features for analyzing protein sequences.

Classes
=======

    * :class:`.MSA` - store MSA data indexed by label

MSA Shuffle
======

    * :func:`.shuffleOMES` - Shuffle the MSA and performe the p-test for OMES
    * :func:`.shuffleMI` - Shuffle the MSA and performe the p-test for MI
    * :func:`.shuffleMIp` - Shuffle the MSA and performe the p-test for MIp
    * :func:`.shuffleAll` - Shuffle the MSA and performe the p-test for all methods

MSA IO
========

    * :func:`.readFasta` - read MSA from FASTA files
    * :func:`.writeFasta` - write MSA to FASTA files

Matrix IO
========

    * :func:`.readMatrix` - read matrix/matrices list from binary file(s)
    * :func:`.writeMatrix` - write matrix to file

Sequence Correlation Analysis
========

    * :func:`.buildMI` - build the MI matrix from MSA
    * :func:`.buildMIp` - build the MIp matrix from MSA
    * :func:`.buildOMES` - build the OMES matrix from MSA
    * :func:`.buildSCA` - build the SCA matrix from MSA
    * :func:`.buildDI` - build the DI(Direct Information) matrix from MSA
    * :func:`.buildDCA` - build the DCA(same as DI) matrix from MSA
    * :func:`.calcMeff`- calc the effective number of sequences
    * :func:`.applyAPC` - apply the APC correction
    * :func:`.applyBND` - apply the BND correction

Sequence
========
    * :func:`.combineMSA` - combine 2 or more MSAs by the labels

"""

__author__ = 'Wenzhi Mao'

__all__ = []


from . import correlation
from .correlation import *
__all__.extend(correlation.__all__)

from . import shuffle
from .shuffle import *
__all__.extend(shuffle.__all__)

from . import msa
from .msa import *
__all__.extend(msa.__all__)

from . import msaio
from .msaio import *
__all__.extend(msaio.__all__)

from . import analysis
from .analysis import *
__all__.extend(analysis.__all__)
