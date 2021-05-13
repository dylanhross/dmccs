#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3
"""

"""
from sqlite3 import connect
from matplotlib import pyplot as plt, rcParams
from numpy import loadtxt, array, linspace
from pickle import load as pload
from scipy.stats import linregress

from DmimData.data import DMD
from metrics import all_metrics


# define a global pRNG seed
SEED = 420

# set the global fontsize on plots
rcParams['font.size'] = 8

# define a query to get structure ids
qry = """
SELECT
    plate_{n}_3d.str_id, pmi1, pmi2, pmi3, rmd02, rmd24, rmd46, rmd68, rmd8p, ccs_avg, annotation, adduct, met_n, mz
FROM
    plate_{n}_md3d
    JOIN plate_{n}_3d
        ON plate_{n}_md3d.str_id = plate_{n}_3d.str_id
    JOIN plate_{n}_id
        ON plate_{n}_3d.ann_id = plate_{n}_id.ann_id
    JOIN plate_{n}
        ON plate_{n}_id.dmim_id = plate_{n}.dmim_id
;"""


def get_indices(cursor, str_id):
    """ fetch the indices of structures with PA/EHS CCS values """
    str_ids = []
    for n in range(1, 8):
        for sid, *_ in cursor.execute(qry.format(n=n)):
            str_ids.append(int(sid))
    idx = []
    for i in range(len(str_ids)):
        if str_ids[i] in str_id:
            idx.append(i)
    return array(idx)


def compare_ccs(y, y_pred, pa, ehs):
    """ make a plot comparing PA/EHS CCS and SVR predicted CCS """

    # perform linear regression on all 3 sets
    svr_slope, svr_intercept, svr_r, svr_p, svr_se = linregress(y_pred, y)
    pa_slope, pa_intercept, pa_r, pa_p, pa_se = linregress(pa, y)
    ehs_slope, ehs_intercept, ehs_r, ehs_p, ehs_se = linregress(ehs, y)

    fig = plt.figure(figsize=(3.3, 3.3))
    ax = fig.add_subplot(111)
    min_y, max_y = min(y), max(y)
    ax.plot([min_y, max_y], [min_y, max_y], 'k-', lw=1, zorder=-1)
    ax.scatter(ehs, y, marker='.', s=2, c='#198C75', edgecolors='none', label='EHS', alpha=1.0)
    ax.scatter(pa, y, marker='.', s=2, c='#0A2F51', edgecolors='none', label='PA', alpha=0.4)
    ax.scatter(y_pred, y, marker='.', s=2, c='#56B870', edgecolors='none', label='SVR', alpha=0.4)
    ehs_s, pa_s, y_pred_s = array(sorted(ehs)), array(sorted(pa)), array(sorted(y_pred))
    ax.plot(ehs_s, ehs_slope * ehs_s + ehs_intercept, ls='--', lw=0.75, c='#198C75')
    ax.plot(pa_s, pa_slope * pa_s + pa_intercept, ls='--', lw=0.75, c='#0A2F51')
    ax.plot(y_pred_s, svr_slope * y_pred_s + svr_intercept, ls='--', lw=0.75, c='#56B870')
    ax.text(0.3, 0.9, r'y = {:.2f}x + {:.2f} $(R^2 = {:.4f})$'.format(ehs_slope, ehs_intercept, ehs_r**2.), 
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=6, 
            c='#198C75', fontweight='bold')
    ax.text(0.3, 0.95, r'y = {:.2f}x + {:.2f} $(R^2 = {:.4f})$'.format(pa_slope, pa_intercept, pa_r**2.), 
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=6, 
            c='#0A2F51', fontweight='bold')
    ax.text(0.29, 0.85, r'y = {:.2f}x + {:.2f} $(R^2 = {:.4f})$'.format(svr_slope, svr_intercept, svr_r**2.), 
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=6, 
            c='#56B870', fontweight='bold')

    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.set_xlabel(r'calc. CCS ($\AA^2$) ', fontsize=8)
    ax.set_ylabel(r'meas. CCS ($\AA^2$)', fontsize=8)
    fig.savefig('PA_EHS_SVR_comparison.png', dpi=400, bbox_inches='tight')
    plt.close()


    # plot the residuals from the linear fits
    svr_res = y - (svr_slope * y_pred + svr_intercept)
    pa_res = y - (pa_slope * pa + pa_intercept)
    ehs_res = y - (ehs_slope * ehs + ehs_intercept)

    fig = plt.figure(figsize=(3.3, 1.3))
    ax = fig.add_subplot(111)
    min_y, max_y = min(y), max(y)
    ax.axhline(0, c='k', ls='--', lw=1, zorder=-1)
    ax.scatter(pa, pa_res, marker='.', s=2, c='#0A2F51', edgecolors='none', label='PA', alpha=0.8)
    ax.scatter(ehs, ehs_res, marker='.', s=2, c='#198C75', edgecolors='none', label='EHS', alpha=0.8)
    ax.scatter(y_pred, svr_res, marker='.', s=2, c='#56B870', edgecolors='none', label='SVR', alpha=0.8)
    ax.set_ylim([-40, 40])
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.set_xlabel(r'calc. CCS ($\AA^2$) ', fontsize=8)
    ax.set_ylabel(r'residual CCS ($\AA^2$)', fontsize=8)
    fig.savefig('PA_EHS_SVR_residuals.png', dpi=400, bbox_inches='tight')
    plt.close()


    # plot the residuals from the linear fits as histograms
    fig = plt.figure(figsize=(3.3, 2.3))
    ax = fig.add_subplot(111)
    bins = linspace(0, 20, 41)
    ax.hist(ehs_res, bins=bins, histtype='stepfilled', color='k', lw=0, alpha=0.1)
    ax.hist(pa_res, bins=bins, histtype='stepfilled', color='k', lw=0, alpha=0.1)
    ax.hist(svr_res, bins=bins, histtype='stepfilled', color='k', lw=0, alpha=0.1)
    ax.hist(ehs_res, bins=bins, histtype='step', color='#198C75', lw=1.5, label='EHS')
    ax.hist(pa_res, bins=bins, histtype='step', color='#0A2F51', lw=1.5, label='PA')
    ax.hist(svr_res, bins=bins, histtype='step', color='#56B870', lw=1.5, label='SVR')
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.set_xlabel(r'residual CCS ($\AA^2$)', fontsize=8)
    ax.set_ylabel('N')
    fig.savefig('PA_EHS_SVR_residuals_hist.png', dpi=400, bbox_inches='tight')
    plt.close()


    




def main():
    """ main execution sequence """
    # set up data interface (custom feature set)
    data = DMD('DMIM_v1.0.db', SEED)
    data.featurize('custom', custom_mqns=['hac', 'c', 'adb', 'asv', 'ctv', 'hbam', 'hbd'], custom_md3ds=['pmi1', 'pmi2', 'pmi3', 'rmd02'])
    data.train_test_split()
    data.center_and_scale()

    # set up a database connection
    con = connect(data.db_path_)
    cur = con.cursor()

    # load the trained SVR
    with open('cust_svr_seed420.pickle', 'rb') as pf:
        svr = pload(pf) 

    # load the PA/EHS CCS
    str_id, pa, ehs = loadtxt('strid_pa_ehs.txt', unpack=True)

    # get the indices of the structures with PA/EHS CCS values
    idx = get_indices(cur, str_id)

    # select out the comparison data
    X_ss = data.SScaler_.transform(data.X_[idx])
    y = data.y_[idx]
    y_pred = svr.predict(X_ss)

    # compare the CCS values
    compare_ccs(y, y_pred, pa, ehs)
    print('SVR predicted')
    #print(all_metrics(y, y_pred))
    print('PA')
    #print(all_metrics(y, pa))
    print('EHS')
    #print(all_metrics(y, ehs))


if __name__ == '__main__':
    main()

