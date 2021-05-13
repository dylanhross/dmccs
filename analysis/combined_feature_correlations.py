#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3
"""

"""

from matplotlib import pyplot as plt, cm, rcParams
from scipy.stats import pearsonr, spearmanr

from DmimData.data import DMD


# define a global pRNG seed
SEED = 420

# set the global fontsize on plots
rcParams['font.size'] = 6


def plot_corr_mats(comb_X_ss):
    """ generate a big figure with a correlation matrix """
    feature_labels = [
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
    ] + ['pmi1', 'pmi2', 'pmi3', 'rmd02', 'rmd24', 'rmd46', 'rmd68', 'rmd8p']
    for a in range(5):
        for b in range(5):
            fig, ax = plt.subplots(10, 10, figsize=(7, 7))
            spec = fig.add_gridspec(ncols=50, nrows=50)
            for i in range(10):
                for j in range(10):
                    ii = a * 10 + i
                    jj = b * 10 + j
                    # Pearson correlation coefficient
                    pr, pv = spearmanr(comb_X_ss.T[ii], comb_X_ss.T[jj])
                    ax[j, i].scatter(comb_X_ss.T[ii], comb_X_ss.T[jj], s=2, marker='.', edgecolors='none', c='b')
                    fw = 'bold' if abs(pr) > 0.8 else 'normal'
                    fs = 'normal' if abs(pr) > 0.8 else 'italic'
                    ax[j, i].text(0.8, 0.8, '{:+.2f}'.format(pr), horizontalalignment='center', 
                                  verticalalignment='center', transform=ax[j, i].transAxes, c='r', fontweight=fw, 
                                  fontstyle=fs)
                    print(pv)
                    if pv < 0.01:
                        ax[j, i].text(0.8, 0.55, '*', horizontalalignment='center', 
                                          verticalalignment='center', transform=ax[j, i].transAxes, c='r', 
                                          fontweight=fw, fontsize=8)
                    for d in ['top', 'right']:
                        ax[i, j].spines[d].set_visible(False)
                    if i == 0:
                        ax[j, i].set_ylabel(feature_labels[jj], fontweight='bold')
                    else:
                        ax[j, i].set_yticks([])
                    if j == 0:
                        ax[j, i].set_xlabel(feature_labels[ii], fontweight='bold')
                        ax[j, i].xaxis.set_label_position('top') 
                    if j != 9:
                        ax[j, i].set_xticks([])
                    

            fig.savefig('corr_mats/DMIM_COMB_corr-mat_{}_{}.png'.format(a, b), dpi=400, bbox_inches='tight')
            plt.close()


def plot_corr_mats2(cust_X_ss):
    """ generate a big figure with a correlation matrix """
    feature_labels = ['hac', 'c', 'adb', 'asv', 'ctv', 'hbam', 'hbd'] + ['pmi1', 'pmi2', 'pmi3', 'rmd02']

    fig, ax = plt.subplots(11, 11, figsize=(7, 7))
    spec = fig.add_gridspec(ncols=11, nrows=11)
    for i in range(11):
        for j in range(11):

            # Pearson correlation coefficient
            pr, pv = spearmanr(cust_X_ss.T[i], cust_X_ss.T[j])
            ax[j, i].scatter(cust_X_ss.T[i], cust_X_ss.T[j], s=1, marker='.', edgecolors='none', c='b')
            fw = 'bold' if abs(pr) > 0.8 else 'normal'
            fs = 'normal' if abs(pr) > 0.8 else 'italic'
            ax[j, i].text(0.8, 0.8, '{:+.2f}'.format(pr), horizontalalignment='center', 
                          verticalalignment='center', transform=ax[j, i].transAxes, c='r', fontweight=fw, 
                          fontstyle=fs)
            print(pv)
            if pv < 0.01:
                ax[j, i].text(0.8, 0.55, '*', horizontalalignment='center', 
                                  verticalalignment='center', transform=ax[j, i].transAxes, c='r', 
                                  fontweight=fw, fontsize=6)
            for d in ['top', 'right']:
                ax[i, j].spines[d].set_visible(False)
            if i == 0:
                ax[j, i].set_ylabel(feature_labels[j], fontweight='bold')
            else:
                ax[j, i].set_yticks([])
            if j == 0:
                ax[j, i].set_xlabel(feature_labels[i], fontweight='bold')
                ax[j, i].xaxis.set_label_position('top') 
            if j != 6:
                ax[j, i].set_xticks([])
                    

    fig.savefig('corr_mats/DMIM_CUST_corr-mat.png', dpi=400, bbox_inches='tight')
    plt.close()



def main():
    """  """ 

    """
    # setup DmimData instance (using MQN + MD3D)
    comb = DMD('DMIM_v1.0.db', SEED)
    comb.featurize('combined')
    comb.train_test_split()
    comb.center_and_scale()
    comb_X_ss = comb.SScaler_.transform(comb.X_)  # scale the MD3D data

    # plot a correlation matrix
    plot_corr_mats(comb_X_ss)
    """

    # setup DmimData instance (using combined MQNs and MD3Ds)
    cust = DMD('DMIM_v1.0.db', SEED)
    cust.featurize('custom', custom_mqns=['hac', 'c', 'adb', 'asv', 'ctv', 'hbam', 'hbd'], custom_md3ds=['pmi1', 'pmi2', 'pmi3', 'rmd02'])
    cust.train_test_split()
    cust.center_and_scale()
    cust_X_ss = cust.SScaler_.transform(cust.X_)

    # plot a correlation matrix
    plot_corr_mats2(cust_X_ss)

if __name__ == '__main__':
    main()
