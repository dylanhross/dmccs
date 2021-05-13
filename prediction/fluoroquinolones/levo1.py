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
    n = int(sys.argv[1])
    
    smis = ['CC1COC3=C2N1C=C(C(=O)C2=CC(=C3N4CC[N+]([H])(CC4)C)F)C(=O)O' for _ in range(n)]
    smis += ['CC1COC3=C2N1C=C(C(=[O+][H])C2=CC(=C3N4CCN(CC4)C)F)C(=O)O' for _ in range(n)]
    smis += ['CC1COC3=C2[N+]1([H])C=C(C(=O)C2=CC(=C3N4CCN(CC4)C)F)C(=O)O' for _ in range(n)]
    smis += ['CC1COC3=C2N1C=C(C(=O)C2=CC(=C3[N+]4([H])CCN(CC4)C)F)C(=O)O' for _ in range(n)]
    structures = []
    for i in range(1, n + 1):
        fname = 'levo_pA_{}.xyzmq'.format(i)
        with open(fname, 'r') as f:
            structures.append(f.read())
    for i in range(1, n + 1):
        fname = 'levo_pB_{}.xyzmq'.format(i)
        with open(fname, 'r') as f:
            structures.append(f.read())
    for i in range(1, n + 1):
        fname = 'levo_pC_{}.xyzmq'.format(i)
        with open(fname, 'r') as f:
            structures.append(f.read())
    for i in range(1, n + 1):
        fname = 'levo_pD_{}.xyzmq'.format(i)
        with open(fname, 'r') as f:
            structures.append(f.read())

    X_cust = featurize(smis, structures, ['hac', 'c', 'adb', 'asv', 'ctv', 'hbam', 'hbd'], ['pmi1', 'pmi2', 'pmi3', 'rmd02'])
    X_mqn = featurize(smis, structures, 'all', [])
    X_md3d = featurize(smis, structures, [], 'all')
    X_comb = featurize(smis, structures, 'all', 'all')

    savetxt('levo_X_CUST.txt', X_cust)
    savetxt('levo_X_MQN.txt', X_mqn)
    savetxt('levo_X_MD3D.txt', X_md3d)
    savetxt('levo_X_COMB.txt', X_comb)
    """
    cust = DMD('DMIM_v1.0.db', SEED)
    cust.featurize('custom', custom_mqns=['hac', 'c', 'asv', 'ctv'], custom_md3ds=['pmi1', 'pmi2', 'pmi3'])
    cust.train_test_split()
    cust.center_and_scale()

    X_ss = cust.SScaler_.transform(X)

    with open('', 'rb') as pf:
        svr = pload(pf)

    print(svr.predict(X_ss))
    """


if __name__ == '__main__':
    main()

