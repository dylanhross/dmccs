#!/usr/local/Cellar/python@3.9/3.9.1_6/bin/python3 
"""

"""

from pickle import load as pload
from numpy import savetxt
import sys

#from DmimData.data import DMD
from helpers import featurize


def main():
    """ main execution sequence """
    n = 4
    fbases = ['q3g', 'q3Pg', 'q4Pg', 'q5g', 'q7g']
    structures, smis = [], []
    for fb in fbases:
        smif = fb + '.smi'
        strucf = 'cactus_' + fb + '.xyzmq'
        with open(smif, 'r') as f:
            smis.append(f.read().strip())
        with open(strucf, 'r') as f:
            structures.append(f.read())

    X = featurize(smis, structures, ['hac', 'c', 'adb', 'asv', 'ctv', 'hbam', 'hbd'], ['pmi1', 'pmi2', 'pmi3', 'rmd02'])
    savetxt('cactus_X.txt', X)


if __name__ == '__main__':
    main()

