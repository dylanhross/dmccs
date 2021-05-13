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
    
    smis = ['CC[N+]1([H])CCN(CC1)C2=C(C=C3C(=C2)N(C=C(C3=O)C(=O)O)C4CC4)F' for _ in range(n)]
    smis += ['CCN1CCN(CC1)C2=C(C=C3C(=C2)N(C=C(C3=[O+][H])C(=O)O)C4CC4)F' for _ in range(n)]
    smis += ['CCN1CCN(CC1)C2=C(C=C3C(=C2)[N+]([H])(C=C(C3=O)C(=O)O)C4CC4)F' for _ in range(n)]
    smis += ['CCN1CC[N+]([H])(CC1)C2=C(C=C3C(=C2)N(C=C(C3=O)C(=O)O)C4CC4)F' for _ in range(n)]
    structures = []
    for i in range(1, n + 1):
        fname = 'enro_pA_{}.xyzmq'.format(i)
        with open(fname, 'r') as f:
            structures.append(f.read())
    for i in range(1, n + 1):
        fname = 'enro_pB_{}.xyzmq'.format(i)
        with open(fname, 'r') as f:
            structures.append(f.read())
    for i in range(1, n + 1):
        fname = 'enro_pC_{}.xyzmq'.format(i)
        with open(fname, 'r') as f:
            structures.append(f.read())
    for i in range(1, n + 1):
        fname = 'enro_pD_{}.xyzmq'.format(i)
        with open(fname, 'r') as f:
            structures.append(f.read())

    X_cust = featurize(smis, structures, ['hac', 'c', 'adb', 'asv', 'ctv', 'hbam', 'hbd'], ['pmi1', 'pmi2', 'pmi3', 'rmd02'])
    X_mqn = featurize(smis, structures, 'all', [])
    X_md3d = featurize(smis, structures, [], 'all')
    X_comb = featurize(smis, structures, 'all', 'all')

    savetxt('enro_X_CUST.txt', X_cust)
    savetxt('enro_X_MQN.txt', X_mqn)
    savetxt('enro_X_MD3D.txt', X_md3d)
    savetxt('enro_X_COMB.txt', X_comb)
    
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

