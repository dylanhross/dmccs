#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3
"""

"""

from numpy import array, abs, argsort, savetxt
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler
from matplotlib import pyplot as plt, cm, rcParams

from DmimData.data import DMD


# define a global pRNG seed
SEED = 420

# set the global fontsize on plots
rcParams['font.size'] = 8


def plot_plsra_c12(c1, c2, c=None, c_min=None, c_max=None, alpha=1., use_ax=None, use_fig=None, s=6, colorbar=False):
    """ sets up PLS-RA projection plots """
    if use_ax is None:
        fig = plt.figure(figsize=(3.5, 3.5))
        ax = fig.add_subplot(111)
        ax.axvline(lw=0.5, c='k', zorder=0)
        ax.axhline(lw=0.5, c='k', zorder=0)
    else: 
        ax = use_ax
        fig = use_fig

    if c_min is None:
        ax.scatter(c1, c2, marker='.', s=s, c=c, edgecolors='none', alpha=alpha)
    else:
        scat = ax.scatter(c1, c2, marker='.', s=s, c=c, cmap=cm.plasma, vmin=c_min, vmax=c_max, edgecolors='none', alpha=alpha)
        if colorbar:
            cax = fig.add_axes([0.2, 0.8, 0.6, 0.05])
            cb = fig.colorbar(scat, cax=cax, orientation='horizontal', ticks=[c_min, c_max])

    if use_ax is None:
        for d in ['top', 'right', 'bottom', 'left']:
            ax.spines[d].set_visible(False)
        ax.set_xlabel('scores[0]', fontsize=8)
        ax.set_ylabel('scores[1]', fontsize=8)
        ax.ticklabel_format(style='sci', scilimits=(0, 0))
        return fig, ax


def compute_mqn_plsra(mqn, ccs):
    """ compute a PLS-RA using the MQN data, return the fitted PLS-RA instance
        make a plot of the projections and save it """
    plsra = PLSRegression(scale=False)  
    projections = plsra.fit_transform(mqn, ccs)
    print('PLS-RA (MQN)')

    # make a plot
    p1, p2 = projections[0].T[0], projections[0].T[1]
    fig, ax = plot_plsra_c12(p1, p2, c=ccs, c_min=120, c_max=240, colorbar=True, s=16)
    fig.savefig('DMIM_MQN_plsra12_cCCS.png', dpi=400, bbox_inches='tight')
    plt.close()
    return plsra


def compute_md3d_plsra(md3d, ccs):
    """ compute a PLS-RA using the MD3D data, return the fitted PCA instance
        make a plot of the projections and save it """
    plsra = PLSRegression(scale=False)  
    projections = plsra.fit_transform(md3d, ccs)
    print('PLS-RA (MD3D)')

    # make a plot
    p1, p2 = projections[0].T[0], projections[0].T[1]
    fig, ax = plot_plsra_c12(p1, p2, c=ccs, c_min=120, c_max=240, colorbar=True, s=8)
    fig.savefig('DMIM_MD3D_plsra12_cCCS.png', dpi=400, bbox_inches='tight')
    plt.close()
    return plsra


def compute_comb_plsra(comb, ccs):
    """ compute a PLS-RA using the MQN + MD3D data, return the fitted PLS-RA instance
        make a plot of the projections and save it """
    plsra = PLSRegression(scale=False)  
    projections = plsra.fit_transform(comb, ccs)
    print('PLS-RA (combined)')

    # make a plot
    p1, p2 = projections[0].T[0], projections[0].T[1]
    fig, ax = plot_plsra_c12(p1, p2, c=ccs, c_min=120, c_max=240, colorbar=True, s=16)
    fig.savefig('DMIM_COMB_plsra12_cCCS.png', dpi=400, bbox_inches='tight')
    plt.close()
    return plsra


def report_mqn_x_loadings(mqn_plsra):
    """ print PCN loadings in order of descending absolute magnitude, make a plot """
    idx = argsort(abs(mqn_plsra.x_loadings_.T[0]))[::-1]
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
    print('feature X loadings')
    print(labels[idx])
    #print(mqn_pca.components_[pc_n - 1][idx])

    # plot the loadings
    fig = plt.figure(figsize=(4, 1.5))
    ax = fig.add_subplot(111)

    y = mqn_plsra.x_loadings_.T[0][idx]
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
    ax.set_ylabel('x loading')

    fig.savefig('DMIM_MQN_plsra_x-loadings.png', dpi=400, bbox_inches='tight')
    plt.close()


def report_md3d_x_loadings(md3d_plsra):
    """ print PCN loadings in order of descending absolute magnitude, make a plot """
    idx = argsort(abs(md3d_plsra.x_loadings_.T[0]))[::-1]
    labels = array(['pmi1', 'pmi2', 'pmi3', 'rmd02', 'rmd24', 'rmd46', 'rmd68', 'rmd8p'])
    print('feature X loadings')
    print(labels[idx])
    #print(md3d_pca.components_[pc_n - 1][idx])

    # plot the loadings
    fig = plt.figure(figsize=(1.5, 1.5))
    ax = fig.add_subplot(111)

    y = md3d_plsra.x_loadings_.T[0][idx]
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
    ax.set_ylabel('x loading')

    fig.savefig('DMIM_MD3D_plsra_x-loadings.png', dpi=400, bbox_inches='tight')
    plt.close()


def report_comb_x_loadings(comb_plsra):
    """ print PCN loadings in order of descending absolute magnitude, make a plot """
    idx = argsort(abs(comb_plsra.x_loadings_.T[0]))[::-1]
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
    print('feature X loadings')
    print(labels[idx])
    #print(mqn_pca.components_[pc_n - 1][idx])

    # plot the loadings
    fig = plt.figure(figsize=(4, 1.5))
    ax = fig.add_subplot(111)

    y = comb_plsra.x_loadings_.T[0][idx]
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
    ax.set_ylabel('x loading')

    fig.savefig('DMIM_COMB_plsra_x-loadings.png', dpi=400, bbox_inches='tight')
    plt.close()

    savetxt('DMIM_COMB_plsra_x-loadings.txt', comb_plsra.x_loadings_.T[0])


def main():
    """ main execution sequence """

    # setup DmimData instance (using MQNs)
    mqn = DMD('DMIM_v1.0.db', SEED)
    mqn.featurize('mqn')
    mqn.train_test_split()
    mqn.center_and_scale()
    mqn_X_ss = mqn.SScaler_.transform(mqn.X_)  # scale the MQN data

    # compute a PCA on the CCSbase MQN data
    mqn_plsra = compute_mqn_plsra(mqn_X_ss, mqn.y_)
    report_mqn_x_loadings(mqn_plsra)

    # setup DmimData instance (using MD3Ds)
    md3d = DMD('DMIM_v1.0.db', SEED)
    md3d.featurize('md3d')
    md3d.train_test_split()
    md3d.center_and_scale()
    md3d_X_ss = md3d.SScaler_.transform(md3d.X_)  # scale the MD3D data

    # compute a PCA on the CCSbase MQN data
    md3d_plsra = compute_md3d_plsra(md3d_X_ss, md3d.y_)
    report_md3d_x_loadings(md3d_plsra)

    # setup DmimData instance (using MQN + MD3D)
    comb = DMD('DMIM_v1.0.db', SEED)
    comb.featurize('combined')
    comb.train_test_split()
    comb.center_and_scale()
    comb_X_ss = comb.SScaler_.transform(comb.X_)  # scale the MD3D data

    # compute a PCA on the CCSbase MQN data
    comb_plsra = compute_comb_plsra(comb_X_ss, comb.y_)
    report_comb_x_loadings(comb_plsra)


if __name__ == '__main__':
    main()
