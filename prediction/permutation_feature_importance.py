#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3
"""

"""
from matplotlib import pyplot as plt, rcParams
from numpy import argsort, array, savetxt
from pickle import load as pload
from sklearn.inspection import permutation_importance

from DmimData.data import DMD


# define a global pRNG seed
SEED = 420

# set the global fontsize on plots
rcParams['font.size'] = 8



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
    ax.set_ylabel('feat. importance')

    fig.savefig('DMIM_COMB_PER_feat_imp.png', dpi=400, bbox_inches='tight')
    plt.close()


def main():
    """ main execution sequence """
    
    data = DMD('DMIM_v1.0.db', SEED)
    data.featurize('combined')
    data.train_test_split()
    data.center_and_scale()
    X = data.SScaler_.transform(data.X_)
    y = data.y_

    with open('comb_svr_seed420.pickle', 'rb') as pf:
        svr = pload(pf) 

    # compute the permutation feature importance
    result = permutation_importance(svr, X, y, n_jobs=-1)
    fi, fi_sd = result.importances_mean, result.importances_std

    # make a plot of the combined feature importances
    report_comb_feat_imp(fi, fi_sd)
    
    # dump the feature importances to file
    savetxt('PER_feature_importances.txt', fi)


if __name__ == '__main__':
    main()