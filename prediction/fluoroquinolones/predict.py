#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3 
"""

"""

from pickle import load as pload
from numpy import loadtxt, savetxt
import sys

from DmimData.data import DMD
#from helpers import featurize


def main():
    """ main execution sequence """
    """
    smis = ['C1CC(=O)NC2=C1C=CC(=C2)OCCCCN3CCN(CC3)C4=C(C(=CC=C4)Cl)Cl', 
            'C1CC(=O)NC2=C1C=CC(=C2)OCCCCN3CCN(CC3)C4=C(C(=CC=C4)Cl)Cl']
    structures = []
    for fname in ['aripip_na_ext3_6-31G.xyzmq', 'aripip_na_fold4_6-31G.xyzmq']:
        with open(fname, 'r'):
            structures.append(f.read())

    X = featurize(smis, structures, ['hac', 'c', 'asv', 'ctv'], ['pmi1', 'pmi2', 'pmi3'])
    """

    cmpd = sys.argv[1]


    X = loadtxt('{}_X_CUST.txt'.format(cmpd))
    cust = DMD('DMIM_v1.0.db', 420)
    cust.featurize('custom', custom_mqns=['hac', 'c', 'adb', 'asv', 'ctv', 'hbam', 'hbd'], custom_md3ds=['pmi1', 'pmi2', 'pmi3', 'rmd02'])
    cust.train_test_split()
    cust.center_and_scale()
    X_ss = cust.SScaler_.transform(X)
    with open('cust_svr_seed420.pickle', 'rb') as pf:
        svr = pload(pf)
    savetxt('{}_y_CUST.txt'.format(cmpd), svr.predict(X_ss))


    X = loadtxt('{}_X_MQN.txt'.format(cmpd))
    mqn = DMD('DMIM_v1.0.db', 420)
    mqn.featurize('mqn')
    mqn.train_test_split()
    mqn.center_and_scale()
    X_ss = mqn.SScaler_.transform(X)
    with open('mqn_svr_seed420.pickle', 'rb') as pf:
        svr = pload(pf)
    savetxt('{}_y_MQN.txt'.format(cmpd), svr.predict(X_ss))


    X = loadtxt('{}_X_MD3D.txt'.format(cmpd))
    md3d = DMD('DMIM_v1.0.db', 420)
    md3d.featurize('md3d')
    md3d.train_test_split()
    md3d.center_and_scale()
    X_ss = md3d.SScaler_.transform(X)
    with open('md3d_svr_seed420.pickle', 'rb') as pf:
        svr = pload(pf)
    savetxt('{}_y_MD3D.txt'.format(cmpd), svr.predict(X_ss))


    X = loadtxt('{}_X_COMB.txt'.format(cmpd))
    comb = DMD('DMIM_v1.0.db', 420)
    comb.featurize('combined')
    comb.train_test_split()
    comb.center_and_scale()
    X_ss = comb.SScaler_.transform(X)
    with open('comb_svr_seed420.pickle', 'rb') as pf:
        svr = pload(pf)
    savetxt('{}_y_COMB.txt'.format(cmpd), svr.predict(X_ss))


if __name__ == '__main__':
    main()

