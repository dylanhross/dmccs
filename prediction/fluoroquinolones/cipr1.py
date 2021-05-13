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
    
    smis = ['C1CC1N2C=C(C(=O)C3=CC(=C(C=C23)N4CC[N+]([H])([H])CC4)F)C(=O)O' for _ in range(n)]
    smis += ['C1CC1N2C=C(C(=[O+][H])C3=CC(=C(C=C23)N4CCNCC4)F)C(=O)O' for _ in range(n)]
    smis += ['C1CC1[N+]2([H])C=C(C(=O)C3=CC(=C(C=C23)N4CCNCC4)F)C(=O)O' for _ in range(n)]
    smis += ['C1CC1N2C=C(C(=O)C3=CC(=C(C=C23)[N+]4([H])CCNCC4)F)C(=O)O' for _ in range(n)]
    structures = []
    for i in range(1, n + 1):
        fname = 'cipr_pA_{}.xyzmq'.format(i)
        with open(fname, 'r') as f:
            structures.append(f.read())
    for i in range(1, n + 1):
        fname = 'cipr_pB_{}.xyzmq'.format(i)
        with open(fname, 'r') as f:
            structures.append(f.read())
    for i in range(1, n + 1):
        fname = 'cipr_pC_{}.xyzmq'.format(i)
        with open(fname, 'r') as f:
            structures.append(f.read())
    for i in range(1, n + 1):
        fname = 'cipr_pD_{}.xyzmq'.format(i)
        with open(fname, 'r') as f:
            structures.append(f.read())

    X_cust = featurize(smis, structures, ['hac', 'c', 'adb', 'asv', 'ctv', 'hbam', 'hbd'], ['pmi1', 'pmi2', 'pmi3', 'rmd02'])
    X_mqn = featurize(smis, structures, 'all', [])
    X_md3d = featurize(smis, structures, [], 'all')
    X_comb = featurize(smis, structures, 'all', 'all')

    savetxt('cipr_X_CUST.txt', X_cust)
    savetxt('cipr_X_MQN.txt', X_mqn)
    savetxt('cipr_X_MD3D.txt', X_md3d)
    savetxt('cipr_X_COMB.txt', X_comb)

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

