#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3 
"""

"""

from pickle import load as pload
from numpy import loadtxt, savetxt

from DmimData.data import DMD
#from helpers import featurize


def main():


    X = loadtxt('X_CUST.txt')
    cust = DMD('DMIM_v1.0.db', 420)
    cust.featurize('custom', custom_mqns=['hac', 'c', 'adb', 'asv', 'ctv', 'hbam', 'hbd'], custom_md3ds=['pmi1', 'pmi2', 'pmi3', 'rmd02'])
    cust.train_test_split()
    cust.center_and_scale()
    X_ss = cust.SScaler_.transform(X)
    with open('cust_svr_seed420.pickle', 'rb') as pf:
        svr = pload(pf)
    savetxt('y_CUST.txt', svr.predict(X_ss))


    X = loadtxt('X_MQN.txt')
    mqn = DMD('DMIM_v1.0.db', 420)
    mqn.featurize('mqn')
    mqn.train_test_split()
    mqn.center_and_scale()
    X_ss = mqn.SScaler_.transform(X)
    with open('mqn_svr_seed420.pickle', 'rb') as pf:
        svr = pload(pf)
    savetxt('y_MQN.txt', svr.predict(X_ss))


    X = loadtxt('X_MD3D.txt')
    md3d = DMD('DMIM_v1.0.db', 420)
    md3d.featurize('md3d')
    md3d.train_test_split()
    md3d.center_and_scale()
    X_ss = md3d.SScaler_.transform(X)
    with open('md3d_svr_seed420.pickle', 'rb') as pf:
        svr = pload(pf)
    savetxt('y_MD3D.txt', svr.predict(X_ss))


    X = loadtxt('X_COMB.txt')
    comb = DMD('DMIM_v1.0.db', 420)
    comb.featurize('combined')
    comb.train_test_split()
    comb.center_and_scale()
    X_ss = comb.SScaler_.transform(X)
    with open('comb_svr_seed420.pickle', 'rb') as pf:
        svr = pload(pf)
    savetxt('y_COMB.txt', svr.predict(X_ss))


if __name__ == '__main__':
    main()

