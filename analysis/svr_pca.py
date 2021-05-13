#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3
"""

"""

from numpy import array, abs, argsort, linspace
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib import pyplot as plt, cm, rcParams
from pickle import load
from scipy.interpolate import griddata

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
    """ compute a PCA using the MQN data, return the fitted PCA instance """
    pca = PCA(n_components=0.95, svd_solver='full', random_state=SEED)  # auto-select the N-components to cover 95% of the variance
    projections = pca.fit_transform(mqn)
    return pca


def compute_md3d_pca(md3d, ccs):
    """ compute a PCA using the MD3D data, return the fitted PCA instance """
    pca = PCA(n_components=0.95, svd_solver='full', random_state=SEED)  # auto-select the N-components to cover 95% of the variance
    projections = pca.fit_transform(md3d)
    return pca


def compute_cust_pca(cust, ccs):
    """ compute a PCA using the custom feature set data, return the fitted PCA instance """
    pca = PCA(n_components=0.95, svd_solver='full', random_state=SEED)  # auto-select the N-components to cover 95% of the variance
    projections = pca.fit_transform(cust)
    return pca


def plot_mqn_svr(mqn, mqn_pca, mqn_svr):
    """ project the MQN SVR support vectors into the PCA space """
    # compute projections of the support vectors
    projections = mqn_pca.transform(mqn.X_train_ss_[mqn_svr.support_])
    proj_all = mqn_pca.transform(mqn.SScaler_.transform(mqn.X_))
    p1, p2 = projections.T[0], projections.T[1]
    p1a, p2a = proj_all.T[0], proj_all.T[1]
    evr1, evr2 = mqn_pca.explained_variance_ratio_[0], mqn_pca.explained_variance_ratio_[1]

    # make a plot (color by CCS)
    fig, ax = plot_pca_pc12(p1a, p2a, evr1=evr1, evr2=evr2, c='k', s=4)
    xi = linspace(-8, 20, 100)
    yi = linspace(-10, 10, 100)
    zi = griddata((p1, p2), mqn.y_train_[mqn_svr.support_], (xi[None,:], yi[:,None]), method='linear')
    ax.contourf(xi, yi, zi, cmap=cm.plasma, alpha=0.5)
    plot_pca_pc12(p1a, p2a, c='k', s=4, use_ax=ax, use_fig=fig)
    plot_pca_pc12(p1, p2, c=mqn.y_train_[mqn_svr.support_], c_min=120, c_max=240, colorbar=True, s=6, use_ax=ax, use_fig=fig)
    fig.savefig('DMIM_MQN-svr_SV_pc12_cCCS.png', dpi=400, bbox_inches='tight')
    plt.close()

    # make a plot (color by COEF)
    c = mqn_svr.dual_coef_.ravel()
    c_min, c_max = min(c), max(c)
    fig, ax = plot_pca_pc12(p1a, p2a, evr1=evr1, evr2=evr2, c='k', s=4)
    xi = linspace(-8, 20, 100)
    yi = linspace(-10, 10, 100)
    zi = griddata((p1, p2), c, (xi[None,:], yi[:,None]), method='linear')
    ax.contourf(xi, yi, zi, cmap=cm.plasma, alpha=0.5)
    plot_pca_pc12(p1a, p2a, c='k', s=4, use_ax=ax, use_fig=fig)
    plot_pca_pc12(p1, p2, c=c, c_min=c_min, c_max=c_max, colorbar=True, s=6, use_ax=ax, use_fig=fig)
    fig.savefig('DMIM_MQN-svr_SV_pc12_cCOEF.png', dpi=400, bbox_inches='tight')
    plt.close()


def plot_md3d_svr(md3d, md3d_pca, md3d_svr):
    """ project the MD3D SVR support vectors into the PCA space """
    # compute projections of the support vectors
    projections = md3d_pca.transform(md3d.X_train_ss_[md3d_svr.support_])
    proj_all = md3d_pca.transform(md3d.SScaler_.transform(md3d.X_))
    p1, p2 = projections.T[0], projections.T[1]
    p1a, p2a = proj_all.T[0], proj_all.T[1]
    evr1, evr2 = md3d_pca.explained_variance_ratio_[0], md3d_pca.explained_variance_ratio_[1]

    # make a plot (color by CCS)
    fig, ax = plot_pca_pc12(p1a, p2a, evr1=evr1, evr2=evr2, c='k', s=4)
    xi = linspace(-8, 20, 100)
    yi = linspace(-10, 10, 100)
    zi = griddata((p1, p2), md3d.y_train_[md3d_svr.support_], (xi[None,:], yi[:,None]), method='linear')
    ax.contourf(xi, yi, zi, cmap=cm.plasma, alpha=0.5)
    plot_pca_pc12(p1a, p2a, c='k', s=4, use_ax=ax, use_fig=fig)
    plot_pca_pc12(p1, p2, c=md3d.y_train_[md3d_svr.support_], c_min=120, c_max=240, colorbar=True, s=6, use_ax=ax, use_fig=fig)
    fig.savefig('DMIM_MD3D-svr_SV_pc12_cCCS.png', dpi=400, bbox_inches='tight')
    plt.close()

    # make a plot (color by COEF)
    c = md3d_svr.dual_coef_.ravel()
    c_min, c_max = min(c), max(c)
    fig, ax = plot_pca_pc12(p1a, p2a, evr1=evr1, evr2=evr2, c='k', s=4)
    xi = linspace(-8, 20, 100)
    yi = linspace(-10, 10, 100)
    zi = griddata((p1, p2), c, (xi[None,:], yi[:,None]), method='linear')
    ax.contourf(xi, yi, zi, cmap=cm.plasma, alpha=0.5)
    plot_pca_pc12(p1a, p2a, c='k', s=4, use_ax=ax, use_fig=fig)
    plot_pca_pc12(p1, p2, c=c, c_min=c_min, c_max=c_max, colorbar=True, s=6, use_ax=ax, use_fig=fig)
    fig.savefig('DMIM_MD3D-svr_SV_pc12_cCOEF.png', dpi=400, bbox_inches='tight')
    plt.close()


def plot_cust_svr(cust, cust_pca, cust_svr):
    """ project the MD3D SVR support vectors into the PCA space """
    # compute projections of the support vectors
    projections = cust_pca.transform(cust.X_train_ss_[cust_svr.support_])
    proj_all = cust_pca.transform(cust.SScaler_.transform(cust.X_))
    p1, p2 = projections.T[0], projections.T[1]
    p1a, p2a = proj_all.T[0], proj_all.T[1]
    evr1, evr2 = cust_pca.explained_variance_ratio_[0], cust_pca.explained_variance_ratio_[1]

    # make a plot (color by CCS)
    fig, ax = plot_pca_pc12(p1a, p2a, evr1=evr1, evr2=evr2, c='k', s=4)
    xi = linspace(-8, 20, 100)
    yi = linspace(-5, 15, 100)
    zi = griddata((p1, p2), cust.y_train_[cust_svr.support_], (xi[None,:], yi[:,None]), method='linear')
    ax.contourf(xi, yi, zi, cmap=cm.plasma, alpha=0.5)
    plot_pca_pc12(p1a, p2a, c='k', s=4, use_ax=ax, use_fig=fig)
    plot_pca_pc12(p1, p2, c=cust.y_train_[cust_svr.support_], c_min=120, c_max=240, colorbar=True, s=6, use_ax=ax, use_fig=fig)
    fig.savefig('DMIM_CUST-svr_SV_pc12_cCCS.png', dpi=400, bbox_inches='tight')
    plt.close()

    # make a plot (color by COEF)
    c = cust_svr.dual_coef_.ravel()
    c_min, c_max = min(c), max(c)
    fig, ax = plot_pca_pc12(p1a, p2a, evr1=evr1, evr2=evr2, c='k', s=4)
    xi = linspace(-8, 20, 100)
    yi = linspace(-5, 15, 100)
    zi = griddata((p1, p2), c, (xi[None,:], yi[:,None]), method='linear')
    ax.contourf(xi, yi, zi, cmap=cm.plasma, alpha=0.5)
    plot_pca_pc12(p1a, p2a, c='k', s=4, use_ax=ax, use_fig=fig)
    plot_pca_pc12(p1, p2, c=c, c_min=c_min, c_max=c_max, colorbar=True, s=6, use_ax=ax, use_fig=fig)
    fig.savefig('DMIM_CUST-svr_SV_pc12_cCOEF.png', dpi=400, bbox_inches='tight')
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

    # project the support vectors into the PCA space
    with open('prediction/mqn_svr_seed420.pickle', 'rb') as pf:
        mqn_svr = load(pf)
    plot_mqn_svr(mqn, mqn_pca, mqn_svr)

    # setup DmimData instance (using MD3Ds)
    md3d = DMD('DMIM_v1.0.db', SEED)
    md3d.featurize('md3d')
    md3d.train_test_split()
    md3d.center_and_scale()
    md3d_X_ss = md3d.SScaler_.transform(md3d.X_)  # scale the MD3D data

    # compute a PCA on the CCSbase MD3D data
    md3d_pca = compute_md3d_pca(md3d_X_ss, md3d.y_)

    # project the support vectors into the PCA space
    with open('prediction/md3d_svr_seed420.pickle', 'rb') as pf:
        md3d_svr = load(pf)
    plot_md3d_svr(md3d, md3d_pca, md3d_svr)

    # setup DmimData instance (using MD3Ds)
    cust = DMD('DMIM_v1.0.db', SEED)
    cust.featurize('custom', custom_mqns=['hac', 'c', 'asv', 'ctv'], custom_md3ds=['pmi1', 'pmi2', 'pmi3'])
    cust.train_test_split()
    cust.center_and_scale()
    cust_X_ss = cust.SScaler_.transform(cust.X_)  # scale the MD3D data

    # compute a PCA on the CCSbase MD3D data
    cust_pca = compute_cust_pca(cust_X_ss, cust.y_)

    # project the support vectors into the PCA space
    with open('prediction/cust_svr_seed420.pickle', 'rb') as pf:
        cust_svr = load(pf)
    plot_cust_svr(cust, cust_pca, cust_svr)



if __name__ == '__main__':
    main()
