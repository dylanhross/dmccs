#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3
"""

"""

from sklearn.svm import SVR
from sklearn.model_selection import GridSearchCV
from pickle import dump

from DmimData.data import DMD
from metrics import all_metrics


# define a global pRNG seed
SEED = 420


def train_svr(X, y):
    """ trains an SVR using X and y data, returns the trained model instance """
    pg = {'C': [10., 100., 1000.], 'gamma': ['scale', 0.001, 0.01, 0.1]}
    gs = GridSearchCV(SVR(cache_size=2048, tol=5e-4, kernel='rbf'), param_grid=pg, n_jobs=-1, cv=5, scoring='neg_mean_squared_error', )
    gs.fit(X, y)
    print('best params:', gs.best_params_)
    return gs.best_estimator_


def main():
    """ main execution sequence """

    from sklearn import __version__ as version

    print('sklearn version:', version)


    # setup DmimData instance (using MQNs)
    mqn = DMD('DMIM_v1.0.db', SEED)
    mqn.featurize('mqn')
    mqn.train_test_split()
    mqn.center_and_scale()
    print('training SVR using MQNs ...')
    mqn_svr = train_svr(mqn.X_train_ss_, mqn.y_train_)
    print('training set:')
    print(all_metrics(mqn.y_train_, mqn_svr.predict(mqn.X_train_ss_)))
    print('test set:')
    print(all_metrics(mqn.y_test_, mqn_svr.predict(mqn.X_test_ss_)))
    print()
    # save the fitted SVR
    with open('mqn_svr_seed420.pickle', 'wb') as pf:
        dump(mqn_svr, pf)
    # and save the scaler
    with open('mqn_scaler_seed420.pickle', 'wb') as pf:
        dump(mqn.SScaler_, pf)


    # setup DmimData instance (using MD3Ds)
    md3d = DMD('DMIM_v1.0.db', SEED)
    md3d.featurize('md3d')
    md3d.train_test_split()
    md3d.center_and_scale()
    print('training SVR using MD3Ds ...')
    md3d_svr = train_svr(md3d.X_train_ss_, md3d.y_train_)
    print('training set:')
    print(all_metrics(md3d.y_train_, md3d_svr.predict(md3d.X_train_ss_)))
    print('test set:')
    print(all_metrics(md3d.y_test_, md3d_svr.predict(md3d.X_test_ss_)))
    print()
    # save the fitted SVR
    with open('md3d_svr_seed420.pickle', 'wb') as pf:
        dump(md3d_svr, pf)
    # and save the scaler
    with open('md3d_scaler_seed420.pickle', 'wb') as pf:
        dump(md3d.SScaler_, pf)

    """
    # setup DmimData instance (using combined MQNs and MD3Ds)
    comb = DMD('DMIM_v1.0.db', SEED)
    comb.featurize('combined')
    comb.train_test_split()
    comb.center_and_scale()
    print('training SVR using combined MQNs and MD3Ds ...')
    comb_svr = train_svr(comb.X_train_ss_, comb.y_train_)
    print('training set:')
    print(all_metrics(comb.y_train_, comb_svr.predict(comb.X_train_ss_)))
    print('test set:')
    print(all_metrics(comb.y_test_, comb_svr.predict(comb.X_test_ss_)))
    print()
    # save the fitted SVR
    with open('comb_svr_seed420.pickle', 'wb') as pf:
        dump(comb_svr, pf)
    """



    # setup DmimData instance (using combined MQNs and MD3Ds)
    cust = DMD('DMIM_v1.0.db', SEED)
    cust.featurize('custom', custom_mqns=['hac', 'c', 'adb', 'asv', 'ctv', 'hbam', 'hbd'], custom_md3ds=['pmi1', 'pmi2', 'pmi3', 'rmd02'])
    cust.train_test_split()
    cust.center_and_scale()
    print('training SVR using custom MQNs and MD3Ds ...')
    cust_svr = train_svr(cust.X_train_ss_, cust.y_train_)
    print('training set:')
    print(all_metrics(cust.y_train_, cust_svr.predict(cust.X_train_ss_)))
    print('test set:')
    print(all_metrics(cust.y_test_, cust_svr.predict(cust.X_test_ss_)))
    print()
    # save the fitted SVR
    with open('cust_svr_seed420.pickle', 'wb') as pf:
        dump(cust_svr, pf)
    # and save the scaler
    with open('cust_scaler_seed420.pickle', 'wb') as pf:
        dump(mqn.SScaler_, pf)
    


if __name__ == '__main__':
    main()

