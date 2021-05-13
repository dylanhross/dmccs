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
        strucf = fb + '.xyzmq'
        with open(smif, 'r') as f:
            smis.append(f.read().strip())
        with open(strucf, 'r') as f:
            structures.append(f.read())

    X_cust = featurize(smis, structures, ['hac', 'c', 'adb', 'asv', 'ctv', 'hbam', 'hbd'], ['pmi1', 'pmi2', 'pmi3', 'rmd02'])
    X_mqn = featurize(smis, structures, 'all', [])
    X_md3d = featurize(smis, structures, [], 'all')
    X_comb = featurize(smis, structures, 'all', 'all')

    savetxt('X_CUST.txt', X_cust)
    savetxt('X_MQN.txt', X_mqn)
    savetxt('X_MD3D.txt', X_md3d)
    savetxt('X_COMB.txt', X_comb)


if __name__ == '__main__':
    main()

