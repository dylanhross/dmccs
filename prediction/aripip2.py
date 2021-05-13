#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3 
"""

"""

from pickle import load as pload
from numpy import loadtxt

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

    X = loadtxt('aripip_X.txt')

    cust = DMD('DMIM_v1.0.db', 420)
    cust.featurize('custom', custom_mqns=['hac', 'c', 'asv', 'ctv'], custom_md3ds=['pmi1', 'pmi2', 'pmi3'])
    cust.train_test_split()
    cust.center_and_scale()

    X_ss = cust.SScaler_.transform(X)

    with open('cust_svr_seed420.pickle', 'rb') as pf:
        svr = pload(pf)

    print(svr.predict(X_ss))


if __name__ == '__main__':
    main()

