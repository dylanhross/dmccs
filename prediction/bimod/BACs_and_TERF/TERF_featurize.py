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
    n = 10
    smis = ['CC(C)(C)C1=CC=C(C=C1)C(CCCN2CCC(CC2)C(C3=CC=CC=C3)(C4=CC=CC=C4)O)O' for _ in range(n)]
    structures = []
    for i in range(1, n + 1):
        fname = 'TERF_c{}.xyzmq'.format(i)
        with open(fname, 'r') as f:
            structures.append(f.read())

    X_cust = featurize(smis, structures, ['hac', 'c', 'adb', 'asv', 'ctv', 'hbam', 'hbd'], ['pmi1', 'pmi2', 'pmi3', 'rmd02'])
    X_mqn = featurize(smis, structures, 'all', [])
    X_md3d = featurize(smis, structures, [], 'all')
    X_comb = featurize(smis, structures, 'all', 'all')

    savetxt('TERF_X_CUST.txt', X_cust)
    savetxt('TERF_X_MQN.txt', X_mqn)
    savetxt('TERF_X_MD3D.txt', X_md3d)
    savetxt('TERF_X_COMB.txt', X_comb)


if __name__ == '__main__':
    main()

