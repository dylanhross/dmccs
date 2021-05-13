#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3
"""

"""

from time import time
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import GridSearchCV
from matplotlib import pyplot as plt, rcParams
from numpy import argsort, array, log10, mean, std, savetxt

from DmimData.data import DMD
from metrics import rmse


# define a global pRNG seed
SEED = 420

# set the global fontsize on plots
rcParams['font.size'] = 8


def trial(X, y, seed):
    """ trains an SVR, returns the RMSE """
    pg = {'n_estimators': [10, 30,100], 'max_depth': [2, 3, 4]}
    est = GradientBoostingRegressor(random_state=seed, tol=1e-2, max_features='auto')
    gs = GridSearchCV(est, param_grid=pg, cv=3, n_jobs=-1, scoring='neg_mean_squared_error')
    gs.fit(X, y)
    #print(gs.best_params_)
    return rmse(y, gs.best_estimator_.predict(X)), gs.best_estimator_.feature_importances_


def report_comb_feat_imp(feat_imp, feat_imp_sd):
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

    # plot the loadings
    fig = plt.figure(figsize=(4, 1.5))
    ax = fig.add_subplot(111)

    y = feat_imp[idx]
    pc = ['b' for _ in y]
    for i in range(len(y)):
        if y[i] < 0.0001:
            pc[i] = 'r'
        
    ax.bar([_ + 1 for _ in range(len(y))], y, color=pc)
    ax.errorbar([_ + 1 for _ in range(len(y))], y, yerr=feat_imp_sd[idx], color='k', fmt='none')
    ax.axhline(0, ls='--', c='k', lw=0.75)
    ax.set_xticks([_ + 1 for _ in range(len(y))])
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.set_xticklabels(labels[idx], rotation='vertical', fontsize=6)
    for xtick, c_ in zip(ax.get_xticklabels(), pc):
        xtick.set_color(c_)
    ax.set_ylabel(r'feat. importance')
    #ax.set_yticks([0, -1, -2, -3, -4])

    fig.savefig('DMIM_COMB_GBR_feat_imp.png', dpi=400, bbox_inches='tight')
    plt.close()


def complete_feature_importance_trials(n_trials):
    """ runs a set of feature importance trials on the complete feature set """
    fi_ = []
    for i in range(n_trials):
        print('trial {}'.format(i + 1))
        t_start = time()
        seed = i + 256
        data = DMD('DMIM_v1.0.db', seed)
        data.featurize('combined')
        data.train_test_split()
        data.center_and_scale()
        rmse, feat_imp = trial(data.X_train_ss_, data.y_train_, seed)
        print('\tRMSE: {:.3f}'.format(rmse))
        print('\t{:.2f} seconds'.format(time() - t_start))
        fi_.append(feat_imp)
    fi_ = array(fi_).T
    fi = array([mean(_) for _ in fi_])
    fi_sd = array([std(_) for _ in fi_])

    # make a plot of the combined feature importances
    report_comb_feat_imp(fi, fi_sd)

    # dump the feature importances to file
    savetxt('GBR_feature_importances.txt', fi)


def main():
    """ main execution sequence """
    
    complete_feature_importance_trials(32)
    

if __name__ == '__main__':
    main()

