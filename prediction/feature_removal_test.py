#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3
"""

"""

from matplotlib import pyplot as plt, rcParams
from numpy import argsort, array, loadtxt, abs
from pickle import load as pload
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVR

from DmimData.data import DMD
from metrics import rmse


# define a global pRNG seed
SEED = 420

# set the global fontsize on plots
rcParams['font.size'] = 8


def train_svr(X, y):
    """ trains an SVR using X and y data, returns the trained model instance """
    pg = {'C': [10., 100., 1000.]}
    gs = GridSearchCV(SVR(cache_size=2048, tol=5e-3, kernel='rbf', gamma='scale'), param_grid=pg, n_jobs=-1, cv=3, scoring='neg_mean_squared_error', )
    gs.fit(X, y)
    return gs.best_estimator_


def separate_features_into_lists(features):
    """ makes separate lists of MQN features and MD3D features """
    mqn = []
    md3d = []
    for feature in features:
        if feature in ['pmi1', 'pmi2', 'pmi3', 'rmd02', 'rmd24', 'rmd46', 'rmd68', 'rmd8p']:
            md3d.append(feature)
        else:
            mqn.append(feature)
    return mqn, md3d


def gen_rmse_data(features):
    """ sequentially removes features in order of increasing feature importance """
    rmses = []
    data = DMD('DMIM_v1.0.db', SEED)
    data.featurize('combined')
    data.train_test_split()
    data.center_and_scale()
    rmses.append(rmse(data.y_train_, train_svr(data.X_train_ss_, data.y_train_).predict(data.X_train_ss_)))

    while len(features) > 1:
        features = features[:-1]
        print(' '.join([_ for _ in features]))
        data = DMD('DMIM_v1.0.db', SEED)
        mqns, md3ds = separate_features_into_lists(features)
        data.featurize('custom', custom_mqns=mqns, custom_md3ds=md3ds)
        data.train_test_split()
        data.center_and_scale()
        rmses.append(rmse(data.y_train_, train_svr(data.X_train_ss_, data.y_train_).predict(data.X_train_ss_)))

    return array(rmses)[::-1]



def report_comb_feat_imp(feat_imp, label):
    """ print feature importances in order of descending absolute magnitude, make a plot """
    idx = argsort(feat_imp)[::-1]
    labels = array([
        'c', 'f', 'cl', 'br',
        'i', 's', 'p', 'an',
        'cn', 'ao', 'co', 'hac',
        'hbam', 'hba', 'hbdm',
        'hbd', 'neg', 'pos',
        'asb', 'adb', 'atb', 'csb',
        'cdb', 'ctb', 'rbc',
        'asv', 'adv', 'atv', 'aqv',
        'cdv', 'ctv', 'cqv', 'r3',
        'r4', 'r5', 'r6', 'r7',
        'r8', 'r9', 'rg10', 
        'afr', 'bfr' 
    ] + ['pmi1', 'pmi2', 'pmi3', 'rmd02', 'rmd24', 'rmd46', 'rmd68', 'rmd8p'])
    #print('feature importance')
    #print(labels[idx])

    # make the RMSE data by sequentially removing features from the end of the list 
    rmse = gen_rmse_data(labels[idx])
    #rmse = [10 for _ in range(10)] + [1 for _ in range(40)]  # dummy RMSE data

    # plot the loadings
    fig = plt.figure(figsize=(4, 1.5))
    ax = fig.add_subplot(111)
    ax2 = ax.twinx()

    y = feat_imp[idx]

    rmse_all = rmse[-1]
    lc = ['b' for _ in range(len(rmse))] + ['k']

    for i in range(len(rmse)):
        if rmse[i] > 1.05 * rmse_all:
            lc[i + 1] = 'r'
    lc[0] = 'r'
        
    
    
    ax.set_xticks([_ for _ in range(len(y) + 1)])
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax2.spines['top'].set_visible(False)

    ax.set_xticklabels(labels[idx].tolist() + ['all'], rotation='vertical', fontsize=6)
    for xtick, c_ in zip(ax.get_xticklabels(), lc + ['k']):
        xtick.set_color(c_)

    ax.set_ylabel(r'RMSE ($\AA^2$)')
    
    ax2.bar([_ for _ in range(len(y))], y, color='grey')
    ax.plot([_ + 1 for _ in range(len(y))], rmse, color='k')
    ax2.set_ylabel('feat. importance', color='grey')
    ax2.tick_params(axis='y', labelcolor='grey')
    ax2.tick_params(axis='y', colors='grey')
    ax2.spines['right'].set_color('grey')

    fig.savefig('DMIM_COMB_{}_seq_feat_rem.png'.format(label), dpi=400, bbox_inches='tight')
    plt.close()


def main():
    """ main execution sequence """
    
    """
    # load the feature importance (GBR)
    print('GBR ...')
    gbr_feat_imp = loadtxt('GBR_feature_importances.txt')
    # make a plot of the combined feature importances
    report_comb_feat_imp(gbr_feat_imp, 'GBR')
    print('done\n')
    """
    

    # load the feature importances (PLSR)
    print('PLSR ... ')
    plsr_feat_imp = abs(loadtxt('DMIM_COMB_plsra_x-loadings.txt'))
    # make a plot of the combined feature importances
    report_comb_feat_imp(plsr_feat_imp, 'PLSR')
    print('done\n')


    # load the feature importance (GBR)
    print('PER ...')
    per_feat_imp = loadtxt('PER_feature_importances.txt')
    # make a plot of the combined feature importances
    report_comb_feat_imp(per_feat_imp, 'PER')
    print('done\n')


if __name__ == '__main__':
    main()

