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
    smis = [
        'CC(C)OC(=O)OC(C)OC(=O)C1=C(CSC2[N+]1([H])C(=O)C2NC(=O)C(=NOC)C3=CSC(=N3)N)COC',
        'CC(C)OC(=O)OC(C)OC(=O)C1=C(CSC2N1C(=O)C2NC(=O)C(=NOC)C3=CSC(=[N+]3[H])N)COC'
    ]
    structures = []
    with open('cefprox_pA.xyzmq', 'r') as f:
        structures.append(f.read())
    with open('cefprox_pB.xyzmq', 'r') as f:
        structures.append(f.read())

    X_cust = featurize(smis, structures, ['hac', 'c', 'adb', 'asv', 'ctv', 'hbam', 'hbd'], ['pmi1', 'pmi2', 'pmi3', 'rmd02'])
    X_mqn = featurize(smis, structures, 'all', [])
    X_md3d = featurize(smis, structures, [], 'all')
    X_comb = featurize(smis, structures, 'all', 'all')

    savetxt('cefprox_X_CUST.txt', X_cust)
    savetxt('cefprox_X_MQN.txt', X_mqn)
    savetxt('cefprox_X_MD3D.txt', X_md3d)
    savetxt('cefprox_X_COMB.txt', X_comb)


if __name__ == '__main__':
    main()

