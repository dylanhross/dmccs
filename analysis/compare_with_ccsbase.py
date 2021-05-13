#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3
"""

"""

from sqlite3 import connect
from numpy import array, arange
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib import pyplot as plt, cm, rcParams
from scipy.optimize import curve_fit

from DmimData.data import DMD


# define a global pRNG seed
SEED = 420

# set the global fontsize on plots
rcParams['font.size'] = 8


def get_ccsbase_mqns():
    """ Fetches all MQN data from CCSbase """
    con = connect('CCSbase_v1.2.db')
    cur = con.cursor()
    qry = """
        SELECT 
            mz, ccs, mqns.* 
        FROM mqns 
            JOIN master 
            ON mqns.g_id = master.g_id
    ;"""
    mz, ccs, mqn = [], [], []
    for mz_, ccs_, gid, *mqn_ in cur.execute(qry).fetchall():
        mz.append(mz_)
        ccs.append(ccs_)
        mqn.append(mqn_)
    return array(mz), array(mqn), array(ccs)


def pf(x, A, B, C):
    """ power function for fitting the CCS vs. m/z data """ 
    return A * (x ** B) + C


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


def plot_mz_ccs(mz, ccs, c, alpha=1., use_ax=None, use_fig=None, s=6):
    """ sets up PCA projection plots """
    if use_ax is None:
        fig = plt.figure(figsize=(3.5, 2.75))
        ax = fig.add_subplot(111)
    else: 
        ax = use_ax
        fig = use_fig

    ax.scatter(mz, ccs, marker='.', s=s, c=c, edgecolors='none', alpha=alpha)
    
    if use_ax is None:
        for d in ['top', 'right']:
            ax.spines[d].set_visible(False)
        #ax.set_xlabel('m/z', fontsize=8)
        ax.set_ylabel(r'CCS ($\AA^2$)', fontsize=8)
        #ax.ticklabel_format(style='sci', scilimits=(0, 0))
        return fig, ax


def ccsbase_pca(ccsbase_mqn_scaled, ccsbase_ccs):
    """ compute a PCA using the MQN data, return the fitted PCA instance
        make a plot of the projections and save it """
    pca = PCA(n_components=0.95, svd_solver='full', random_state=SEED)  # auto-select the N-components to cover 95% of the variance
    projections = pca.fit_transform(ccsbase_mqn_scaled)
    print('PCA (CCSbase)')
    print(pca.n_components_, 'components required to explain 95% of the variance')

    # make a plot
    p1, p2 = projections.T[0], projections.T[1]
    evr1, evr2 = pca.explained_variance_ratio_[0], pca.explained_variance_ratio_[1]
    fig, ax = plot_pca_pc12(p1, p2, evr1=evr1, evr2=evr2, c=ccsbase_ccs, c_min=100, c_max=340, colorbar=True)
    fig.savefig('CCSbase_pc12_cCCS.png', dpi=400, bbox_inches='tight')
    plt.close()
    return pca


def plot_dmim_on_ccsbase_pca(dmdata, ccsbase_mz, ccsbase_ccs, ccsbase_mqn, ccsbase_ss, ccsbase_pca):
    """ plot projections of DMIM data using PCA computed on CCSbase, color by CCS, parents/metabolites """
    # get and transform all of the data
    ccsb_proj = ccsbase_pca.transform(ccsbase_ss.transform(ccsbase_mqn))
    ccsb_p1, ccsb_p2 = ccsb_proj.T[0], ccsb_proj.T[1]
    evr1, evr2 = ccsbase_pca.explained_variance_ratio_[0], ccsbase_pca.explained_variance_ratio_[1]
    dmd_proj = ccsbase_pca.transform(ccsbase_ss.transform(dmdata.X_))
    dmd_p1, dmd_p2 = dmd_proj.T[0], dmd_proj.T[1]

    # generate and save the plot (color by CCS)
    fig, ax = plot_pca_pc12(ccsb_p1, ccsb_p2, evr1=evr1, evr2=evr2, c='grey', s=4)
    plot_pca_pc12(dmd_p1, dmd_p2, c=dmdata.y_, c_min=120, c_max=240, use_ax=ax, use_fig=fig, s=8, colorbar=True)
    fig.savefig('DMIM_CCSbase_pc12_cCCS.png', dpi=400, bbox_inches='tight')
    plt.close()

    # generate and save the plot (color by parents/metabolites)
    fig, ax = plot_pca_pc12(ccsb_p1, ccsb_p2, evr1=evr1, evr2=evr2, c='grey', s=4)
    c_pm = ['b' if _ == 0 else 'r' for _ in dmdata.met_n_]
    plot_pca_pc12(dmd_p1, dmd_p2, c=c_pm, use_ax=ax, use_fig=fig, s=8)
    fig.savefig('DMIM_CCSbase_pc12_cMET.png', dpi=400, bbox_inches='tight')
    plt.close()

    # generate and save a plot of CCS vs. m/z of the DMIM data overlaid on the CCSbase data
    fig, ax = plot_mz_ccs(ccsbase_mz, ccsbase_ccs, 'grey', s=12, alpha=0.2)
    par = array([(x, y) for x, y, z in zip(dmdata.mz_, dmdata.y_, dmdata.met_n_) if z == 0]).T
    met = array([(x, y) for x, y, z in zip(dmdata.mz_, dmdata.y_, dmdata.met_n_) if z > 0]).T
    plot_mz_ccs(*par, '#0000FF', s=5, use_ax=ax, use_fig=fig)
    plot_mz_ccs(*met, '#FF0000', s=5, use_ax=ax, use_fig=fig)
    ax.set_xlim([0, 1000])
    ax.set_ylim([100, 350])
    # add fitted curves to the plot
    par_opt, _ = curve_fit(pf, *par, p0=(1., 1., 0.))
    met_opt, _ = curve_fit(pf, *met, p0=(1., 1., 0.))
    x = arange(100, 900, 0.1)
    ax.plot(x, pf(x, *par_opt), '--', c='chartreuse', lw=1)
    ax.plot(x, pf(x, *met_opt), '--', c='darkorange', lw=1)
    fig.savefig('DMIM_CCSbase_mz-ccs_cMET.png', dpi=400, bbox_inches='tight')
    plt.close()

    # generate a plot of residual CCS vs. m/z for just the fitted parent/metabolite data
    fig = plt.figure(figsize=(3.5, 0.75))
    ax = fig.add_subplot(111)
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.set_xlabel('m/z', fontsize=8)
    #ax.set_ylabel(r'resid. CCS ($\AA^2$)', fontsize=8)
    par_resid = pf(par[0], *par_opt) - par[1]
    met_resid = pf(met[0], *met_opt) - met[1]
    ax.scatter(par[0], par_resid, marker='.', s=5, c='#0000FF', edgecolors='none')
    ax.scatter(met[0], met_resid, marker='.', s=5, c='#FF0000', edgecolors='none')
    ax.axhline(0, c='k', ls='--', lw=0.75, zorder=-1)
    ax.set_xlim([0, 1000])
    ax.set_ylim([-40, 40])
    # add fitted curves to the plot
    fig.savefig('DMIM_CCSbase_mz-ccs_resid_cMET.png', dpi=400, bbox_inches='tight')
    plt.close()



def main():
    """ main execution sequence """

    # get MQN data from CCSbase, make a scaled version
    ccsb_mz, ccsb_mqn, ccsb_ccs = get_ccsbase_mqns()
    ccsb_ss = StandardScaler()
    ccsb_mqn_ss = ccsb_ss.fit_transform(ccsb_mqn)

    # compute a PCA on the CCSbase MQN data
    ccsb_pca = ccsbase_pca(ccsb_mqn_ss, ccsb_ccs)

    # load DMIM data and compute projections with the CCSbase PCA, make a plot
    dmd = DMD('DMIM_v1.0.db', SEED)
    dmd.featurize('mqn')
    
    # plot the projections of the DMIM data with the PCA computed on CCSbase
    plot_dmim_on_ccsbase_pca(dmd, ccsb_mz, ccsb_ccs, ccsb_mqn, ccsb_ss, ccsb_pca)




if __name__ == '__main__':
    main()
