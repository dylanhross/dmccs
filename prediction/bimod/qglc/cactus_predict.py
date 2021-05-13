#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3 
"""

"""

from pickle import load as pload
from numpy import loadtxt, savetxt

from DmimData.data import DMD
#from helpers import featurize


def main():


    X = loadtxt('cactus_X.txt')

    cust = DMD('DMIM_v1.0.db', 420)
    cust.featurize('custom', custom_mqns=['hac', 'c', 'adb', 'asv', 'ctv', 'hbam', 'hbd'], custom_md3ds=['pmi1', 'pmi2', 'pmi3', 'rmd02'])
    cust.train_test_split()
    cust.center_and_scale()

    X_ss = cust.SScaler_.transform(X)

    with open('cust_svr_seed420.pickle', 'rb') as pf:
        svr = pload(pf)

    savetxt('cactus_y.txt', svr.predict(X_ss))


if __name__ == '__main__':
    main()

