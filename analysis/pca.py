#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3
"""

"""

from numpy import array, abs, argsort
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib import pyplot as plt, cm, rcParams

from DmimData.data import DMD


# define a global pRNG seed
SEED = 420

# set the global fontsize on plots
rcParams['font.size'] = 8


def plot_pca_pc12(pc1, pc2, evr1=None, evr2=None, c=None, c_min=None, c_max=None, alpha=1., use_ax=None, use_fig=None, s=6, colorbar=False):
    """ sets up PCA projection plots """
    if use_ax is None:
        fig = plt.figure(figsize=(3.5, 3.5))
        ax = fig.add_subplot(111)
        ax.axvline(lw=0.5, c='k', zorder=0)
        ax.axhline(lw=0.5, c='k', zorder=0)
    else: 
        ax = use_ax
        fig = use_fig

    if c_min is None:
        ax.scatter(pc1, pc2, marker='.', s=s, c=c, edgecolors='none', alpha=alpha)
    else:
        scat = ax.scatter(pc1, pc2, marker='.', s=s, c=c, cmap=cm.plasma, vmin=c_min, vmax=c_max, edgecolors='none', alpha=alpha)
        if colorbar:
            cax = fig.add_axes([0.2, 0.8, 0.6, 0.05])
            cb = fig.colorbar(scat, cax=cax, orientation='horizontal', ticks=[c_min, c_max])

    if use_ax is None:
        for d in ['top', 'right', 'bottom', 'left']:
            ax.spines[d].set_visible(False)
        ax.set_xlabel('PC1 ({:.1f} %)'.format(100. * evr1), fontsize=8)
        ax.set_ylabel('PC2 ({:.1f} %)'.format(100. * evr2), fontsize=8)
        ax.ticklabel_format(style='sci', scilimits=(0, 0))
        return fig, ax


def compute_mqn_pca(mqn, ccs):
    """ compute a PCA using the MQN data, return the fitted PCA instance
        make a plot of the projections and save it """
    pca = PCA(n_components=0.95, svd_solver='full', random_state=SEED)  # auto-select the N-components to cover 95% of the variance
    projections = pca.fit_transform(mqn)
    print('PCA (MQN)')
    print(pca.n_components_, 'components required to explain 95% of the variance')

    # make a plot
    p1, p2 = projections.T[0], projections.T[1]
    evr1, evr2 = pca.explained_variance_ratio_[0], pca.explained_variance_ratio_[1]
    fig, ax = plot_pca_pc12(p1, p2, evr1=evr1, evr2=evr2, c=ccs, c_min=120, c_max=240, colorbar=True, s=16)
    fig.savefig('DMIM_MQN_pc12_cCCS.png', dpi=400, bbox_inches='tight')
    plt.close()
    return pca


def compute_md3d_pca(md3d, ccs):
    """ compute a PCA using the MD3D data, return the fitted PCA instance
        make a plot of the projections and save it """
    pca = PCA(n_components=0.95, svd_solver='full', random_state=SEED)  # auto-select the N-components to cover 95% of the variance
    projections = pca.fit_transform(md3d)
    print('PCA (MD3D)')
    print(pca.n_components_, 'components required to explain 95% of the variance')

    # make a plot
    p1, p2 = projections.T[0], projections.T[1]
    evr1, evr2 = pca.explained_variance_ratio_[0], pca.explained_variance_ratio_[1]
    fig, ax = plot_pca_pc12(p1, p2, evr1=evr1, evr2=evr2, c=ccs, c_min=120, c_max=240, colorbar=True, s=8)
    fig.savefig('DMIM_MD3D_pc12_cCCS.png', dpi=400, bbox_inches='tight')
    plt.close()
    return pca


def compute_comb_pca(comb, ccs):
    """ compute a PCA using the MQN + MD3D data, return the fitted PCA instance
        make a plot of the projections and save it """
    pca = PCA(n_components=0.95, svd_solver='full', random_state=SEED)  # auto-select the N-components to cover 95% of the variance
    projections = pca.fit_transform(comb)
    print('PCA (combined)')
    print(pca.n_components_, 'components required to explain 95% of the variance')

    # make a plot
    p1, p2 = projections.T[0], projections.T[1]
    evr1, evr2 = pca.explained_variance_ratio_[0], pca.explained_variance_ratio_[1]
    fig, ax = plot_pca_pc12(p1, p2, evr1=evr1, evr2=evr2, c=ccs, c_min=120, c_max=240, colorbar=True, s=8)
    fig.savefig('DMIM_COMB_pc12_cCCS.png', dpi=400, bbox_inches='tight')
    plt.close()
    return pca


def compute_cust_pca(cust, ccs):
    """ compute a PCA using the customized MQN + MD3D data, return the fitted PCA instance
        make a plot of the projections and save it """
    pca = PCA(n_components=0.95, svd_solver='full', random_state=SEED)  # auto-select the N-components to cover 95% of the variance
    projections = pca.fit_transform(cust)
    print('PCA (custom)')
    print(pca.n_components_, 'components required to explain 95% of the variance')

    # make a plot
    p1, p2 = projections.T[0], projections.T[1]
    evr1, evr2 = pca.explained_variance_ratio_[0], pca.explained_variance_ratio_[1]
    fig, ax = plot_pca_pc12(p1, p2, evr1=evr1, evr2=evr2, c=ccs, c_min=120, c_max=240, colorbar=True, s=8)
    fig.savefig('DMIM_CUST_pc12_cCCS.png', dpi=400, bbox_inches='tight')
    plt.close()
    return pca


def report_mqn_pcN_loadings(mqn_pca, pc_n):
    """ print PCN loadings in order of descending absolute magnitude, make a plot """
    idx = argsort(abs(mqn_pca.components_[pc_n - 1]))[::-1]
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
    ])
    print('feature loadings (PC{})'.format(pc_n))
    print(labels[idx])
    #print(mqn_pca.components_[pc_n - 1][idx])

    # plot the loadings
    fig = plt.figure(figsize=(4, 1.5))
    ax = fig.add_subplot(111)

    y = mqn_pca.components_[pc_n - 1][idx]
    pc = ['b' for _ in y]
    for i in range(len(y)):
        if y[i] < 0:
            pc[i] = 'r'
        
    ax.bar([_ + 1 for _ in range(len(y))], y, color=pc)
    ax.axhline(0, ls='--', c='k', lw=0.75)
    ax.set_xticks([_ + 1 for _ in range(len(y))])
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.set_xticklabels(labels[idx], rotation='vertical', fontsize=6)
    for xtick, c_ in zip(ax.get_xticklabels(), pc):
        xtick.set_color(c_)
    ax.set_ylabel('PC{} loading'.format(pc_n))

    fig.savefig('DMIM_MQN_pc{}-loadings.png'.format(pc_n), dpi=400, bbox_inches='tight')
    plt.close()


def report_md3d_pcN_loadings(md3d_pca, pc_n):
    """ print PCN loadings in order of descending absolute magnitude, make a plot """
    idx = argsort(abs(md3d_pca.components_[pc_n - 1]))[::-1]
    labels = array(['pmi1', 'pmi2', 'pmi3', 'rmd02', 'rmd24', 'rmd46', 'rmd68', 'rmd8p'])
    print('feature loadings (PC{})'.format(pc_n))
    print(labels[idx])
    #print(md3d_pca.components_[pc_n - 1][idx])

    # plot the loadings
    fig = plt.figure(figsize=(1.5, 1.5))
    ax = fig.add_subplot(111)

    y = md3d_pca.components_[pc_n - 1][idx]
    pc = ['b' for _ in y]
    for i in range(len(y)):
        if y[i] < 0:
            pc[i] = 'r'
        
    ax.bar([_ + 1 for _ in range(len(y))], y, color=pc)
    ax.axhline(0, ls='--', c='k', lw=0.75)
    ax.set_xticks([_ + 1 for _ in range(len(y))])
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.set_xticklabels(labels[idx], rotation='vertical', fontsize=6)
    for xtick, c_ in zip(ax.get_xticklabels(), pc):
        xtick.set_color(c_)
    ax.set_ylabel('PC{} loading'.format(pc_n))

    fig.savefig('DMIM_MD3D_pc{}-loadings.png'.format(pc_n), dpi=400, bbox_inches='tight')
    plt.close()


def report_comb_pcN_loadings(comb_pca, pc_n):
    """ print PCN loadings in order of descending absolute magnitude, make a plot """
    idx = argsort(abs(comb_pca.components_[pc_n - 1]))[::-1]
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
    print('feature loadings (PC{})'.format(pc_n))
    print(labels[idx])
    #print(mqn_pca.components_[pc_n - 1][idx])

    # plot the loadings
    fig = plt.figure(figsize=(4, 1.5))
    ax = fig.add_subplot(111)

    y = comb_pca.components_[pc_n - 1][idx]
    pc = ['b' for _ in y]
    for i in range(len(y)):
        if y[i] < 0:
            pc[i] = 'r'
        
    ax.bar([_ + 1 for _ in range(len(y))], y, color=pc)
    ax.axhline(0, ls='--', c='k', lw=0.75)
    ax.set_xticks([_ + 1 for _ in range(len(y))])
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.set_xticklabels(labels[idx], rotation='vertical', fontsize=6)
    for xtick, c_ in zip(ax.get_xticklabels(), pc):
        xtick.set_color(c_)
    ax.set_ylabel('PC{} loading'.format(pc_n))

    fig.savefig('DMIM_COMB_pc{}-loadings.png'.format(pc_n), dpi=400, bbox_inches='tight')
    plt.close()


def main():
    """ main execution sequence """

    # setup DmimData instance (using MQNs)
    mqn = DMD('DMIM_v1.0.db', SEED)
    mqn.featurize('mqn')
    mqn.train_test_split()
    mqn.center_and_scale()
    mqn_X_ss = mqn.SScaler_.transform(mqn.X_)  # scale the MQN data

    # compute a PCA on the CCSbase MQN data
    mqn_pca = compute_mqn_pca(mqn_X_ss, mqn.y_)
    report_mqn_pcN_loadings(mqn_pca, 1)
    report_mqn_pcN_loadings(mqn_pca, 2)
    report_mqn_pcN_loadings(mqn_pca, 3)

    # setup DmimData instance (using MD3Ds)
    md3d = DMD('DMIM_v1.0.db', SEED)
    md3d.featurize('md3d')
    md3d.train_test_split()
    md3d.center_and_scale()
    md3d_X_ss = md3d.SScaler_.transform(md3d.X_)  # scale the MD3D data

    # compute a PCA on the CCSbase MQN data
    md3d_pca = compute_md3d_pca(md3d_X_ss, md3d.y_)
    report_md3d_pcN_loadings(md3d_pca, 1)
    report_md3d_pcN_loadings(md3d_pca, 2)

    # setup DmimData instance (using MQN + MD3D)
    comb = DMD('DMIM_v1.0.db', SEED)
    comb.featurize('combined')
    comb.train_test_split()
    comb.center_and_scale()
    comb_X_ss = comb.SScaler_.transform(comb.X_)  # scale the MD3D data

    # compute a PCA on the CCSbase MQN data
    comb_pca = compute_comb_pca(comb_X_ss, comb.y_)
    report_comb_pcN_loadings(comb_pca, 1)

    # setup DmimData instance (using MQN + MD3D)
    cust = DMD('DMIM_v1.0.db', SEED)
    cust.featurize('custom', custom_mqns=['hac', 'c', 'asv', 'ctv'], custom_md3ds=['pmi1', 'pmi2', 'pmi3'])
    cust.train_test_split()
    cust.center_and_scale()
    cust_X_ss = cust.SScaler_.transform(cust.X_)  # scale the MD3D data

    # compute a PCA on the CCSbase MQN data
    cust_pca = compute_cust_pca(cust_X_ss, cust.y_)


if __name__ == '__main__':
    main()
