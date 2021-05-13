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
    smis = [
        'C[N+](C)(C)CCCC',
        'C[N+](C)(C)CCCCCC',
        'C[N+](C)(C)CCCCCCCC',
        'C[N+](C)(C)CCCCCCCCCC',
        'C[N+](C)(C)CCCCCCCCCCCC',
        'C[N+](C)(C)CCCCCCCCCCCCCC',
        'C[N+](C)(C)CCCCCCCCCCCCCCCC'
    ]
    structures = []
    for i in range(4, 17, 2):
        fname = 'C{:02d}_ext.xyzmq'.format(i)
        with open(fname, 'r') as f:
            structures.append(f.read())

    X_cust = featurize(smis, structures, ['hac', 'c', 'adb', 'asv', 'ctv', 'hbam', 'hbd'], ['pmi1', 'pmi2', 'pmi3', 'rmd02'])
    X_mqn = featurize(smis, structures, 'all', [])
    X_md3d = featurize(smis, structures, [], 'all')
    X_comb = featurize(smis, structures, 'all', 'all')

    savetxt('extended_X_CUST.txt', X_cust)
    savetxt('extended_X_MQN.txt', X_mqn)
    savetxt('extended_X_MD3D.txt', X_md3d)
    savetxt('extended_X_COMB.txt', X_comb)



if __name__ == '__main__':
    main()

