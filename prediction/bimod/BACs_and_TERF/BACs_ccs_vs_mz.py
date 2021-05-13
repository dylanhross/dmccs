#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3
"""
"""


from matplotlib import pyplot as plt
import json
import numpy as np
from scipy.optimize import curve_fit
import sys


def compaction_factor(parent_mass, parent_ccs, metabolite_mass, metabolite_ccs):
    ccs_ratio = parent_ccs  / metabolite_ccs
    mass_ratio = (parent_mass / metabolite_mass)**(2. / 3.)
    return ccs_ratio / mass_ratio


def main():

    # set up plot fonts
    from matplotlib import rcParams
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Helvetica', 'Arial']
    rcParams['font.size'] = 8



    fset = sys.argv[1]

    """
        CCS vs. m/z plot of measured and TJM CCS
    """
    # set up a power fit on the big drugs dataset
    with open('hine0817.json', 'r') as j:
        mz, ccs = np.array([[entry['mz'], entry['ccs']] for entry in json.load(j)]).T

    def pfunc(x, A, B, C):
        return A * (x + B)**C

    opt, cv = curve_fit(pfunc, mz, ccs, p0=(1., 0., 0.5))

    x = np.linspace(min(mz), max(mz), 500)


    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111)

    # plot the power fit on the drugs dataset in the background
    ax.fill_between(x, 0.9 * pfunc(x, *opt), 1.1 * pfunc(x, *opt), facecolor='c', alpha=0.15)
    ax.plot(x, 0.9 * pfunc(x, *opt), c='c', ls='--', lw=1, alpha=0.3)
    ax.plot(x, 1.1 * pfunc(x, *opt), c='c', ls='--', lw=1, alpha=0.3)
    ax.plot(x, pfunc(x, *opt), c='c', ls='-', lw=2, alpha=0.3)


    # observed
    """
    bacs = [[192.1765, 220.2098, 248.2418, 276.2810, 304.3118, 332.3402, 360.3755],  
            [146.03,   158.20,   171.42,   183.78,   193.77,   203.07,   212.31  ]]
    bacs_O = [[208.1724, 236.2084, 264.2374, 292.2759, 320.3052, 348.3354, 376.3714],
              [147.80,   159.51,   170.87,   179.55,   187.83,   195.86,   203.62  ]]
    """
    # use average values from 3 replicates
    bacs_3reps = [
        [192.1838,  144.96,  0.76],
        [208.1795,  148.19,  0.33],
        [220.2167,  157.07,  0.80],
        [236.2128,  158.92,  0.42],
        [248.2505,  170.09,  0.94],
        [264.2128,  170.14,  0.52],
        [276.2869,  182.24,  1.09],
        [292.2815,  178.99,  0.42],
        [304.3187,  193.19,  0.41],
        [320.3119,  187.55,  0.37],
        [332.3427,  203.67,  0.54],
        [348.3444,  196.33,  0.48],
        [360.3843,  212.76,  0.34],
        [376.3800,  205.00,  1.08]
    ]
    bacs, bacs_O, bacs_err, bacs_O_err  = [[], []], [[], []], [], []
    for i in range(len(bacs_3reps)):
        mz, ccs, e = bacs_3reps[i]
        if i % 2 == 0:
            # parent
            bacs[0].append(mz)
            bacs[1].append(ccs)
            bacs_err.append(e)
        else:
            # metabolite
            bacs_O[0].append(mz)
            bacs_O[1].append(ccs)
            bacs_O_err.append(e)

    # compute weighted average of predicted CCS
    wavg, wavgO = [], []
    for i in range(4, 17, 2):
        a = np.loadtxt('C{:02d}_y_{}.txt'.format(i, fset))
        b = np.loadtxt('C{:02d}_cluster_sizewt.txt'.format(i), usecols=(1))
        c = np.loadtxt('C{:02d}OH_y_{}.txt'.format(i, fset))
        d = np.loadtxt('C{:02d}OH_cluster_sizewt.txt'.format(i), usecols=(1))
        wavg.append(sum(a * b) / sum(b))
        wavgO.append(sum(c * d) / sum(d))

    print(wavg)
    print(wavgO)


    # theoretical
    tjm = [[192.1765, 220.2098, 248.2418, 276.2810, 304.3118, 332.3402, 360.3755],  
           wavg]
    tjm_O = [[208.1724, 236.2084, 264.2374, 292.2759, 320.3052, 348.3354, 376.3714],
             wavgO]
    maxs = [[192.1765, 220.2098, 248.2418, 276.2810, 304.3118, 332.3402, 360.3755],  
            np.loadtxt('extended_y_{}.txt'.format(fset))]


    # compute compaction factors for all of the metabolites
    meas_cf = []
    for parent_mass, parent_ccs, metabolite_mass, metabolite_ccs in zip(*bacs, *bacs_O):
        meas_cf.append(compaction_factor(parent_mass, parent_ccs, metabolite_mass, metabolite_ccs))

    pred_cf = []
    for parent_mass, parent_ccs, metabolite_mass, metabolite_ccs in zip(*tjm, *tjm_O):
        pred_cf.append(compaction_factor(parent_mass, parent_ccs, metabolite_mass, metabolite_ccs))

    meas_cf = ['{:.3f}'.format(_) for _ in meas_cf]
    pred_cf = ['{:.3f}'.format(_) for _ in pred_cf]

    print(meas_cf)
    print(pred_cf)

    ms, mew = 6, 1.25
    ax.plot(*bacs, 'ko', ms=ms, label='parents (exp.)', fillstyle='none', mew=mew)
    #ax.errorbar(*bacs, yerr=bacs_err, linestyle='none', c='k')
    ax.plot(*bacs_O, 'k^', ms=ms, label='+O metab. (exp.)', fillstyle='none', mew=mew)
    for x, y, s in zip(*bacs_O, meas_cf):
        ax.text(x + 5, y, s, c='k', fontsize=7)
    #ax.errorbar(*bacs_O, yerr=bacs_O_err, linestyle='none', c='k')
    ax.plot(*tjm, 'bo', ms=ms, label='parents (pred.)', mew=mew)
    #ax.errorbar(*tjm, yerr=tjm_err, linestyle='none', c='b')
    ax.plot(*tjm_O, 'b^', ms=ms, label='+O metab.(pred.)', mew=mew)
    i = 0
    for x, y, s in zip(*tjm_O, pred_cf):
        if fset in ['MQN', 'CUST']:
            yoff = -7 if i == 0 else -5 if i == 0 else -3
        elif fset == 'MD3D':
            yoff = -4
        elif fset == 'COMB':
            yoff = 7 if i <= 3 else -5
        i += 1
        ax.text(x + 5, y + yoff, s, c='b', fontsize=7)
    #ax.errorbar(*tjm_O, yerr=tjm_O_err, linestyle='none', c='b')
    ax.plot(*maxs, 'gx', ms=ms, label='parents \n(extended, pred.)', mew=mew, fillstyle='none')

    ax.set_xlim([150, 400])
    ax.set_ylim([115, 235])


    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.set_ylabel(r'CCS ($\AA^2$)')
    ax.set_xlabel('m/z')
    ax.legend(fontsize=8, loc='upper left', bbox_to_anchor=(0.05, 0.98))
    plt.tight_layout()

    # save figure
    png = 'BACs_{}_ccsvsmz_comp.png'.format(fset)
    plt.savefig(png, dpi=400, bbox_inches='tight')

    #plt.show()
    plt.close()



    """
            Ratio of parent to metabolite CCS, measured and calculated
    """

    """
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111)

    # compute ratios
    n = np.array([4, 6, 8, 10, 12, 14, 16])
    meas_ratio = [a / b for a, b in zip(bacs[1], bacs_O[1])]
    tjm_ratio = [a / b for a, b in zip(tjm[1], tjm_O[1])]
    tjm_ratio_err = [p / m * np.sqrt((pe / p)**2. + (me / m)**2.) for p, m, pe, me in zip(tjm[1],tjm_O[1], tjm_err, tjm_O_err)]

    def lin(x, m, b):
        return m * x + b
    mopt, mcov = curve_fit(lin, n, meas_ratio)
    topt, tcov = curve_fit(lin, n[:-1], tjm_ratio[:-1])

    ###

    ms, mew = 7, 1.25
    ax.axhline(1., ls='-', lw=1.5, c='grey')
    ax.plot(n, meas_ratio, 'ko', ms=ms, label='exp.', fillstyle='none', mew=mew)
    ax.plot(n, lin(n, *mopt), 'k--', linewidth=1.25)
    ax.plot(n, tjm_ratio, 'ko', ms=ms, label='theo.', mew=mew)
    ax.plot(n[:-1], lin(n[:-1], *topt), 'k-', linewidth=1.25)

    #ax.errorbar(n, tjm_ratio, yerr=tjm_ratio_err, linestyle='none', c='k')


    #ax.set_xlim([150, 400])
    ax.set_ylim([0.98, 1.1])

    ax.set_xticks(n)


    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.set_ylabel('parent CCS : +O metab. CCS')
    ax.set_xlabel('BAC chain length')
    ax.legend(fontsize=10)
    plt.tight_layout()

    # save figure
    png = 'BACs_ratio_comp.png'
    plt.savefig(png, dpi=400, bbox_inches='tight')

    #plt.show()
    plt.close()
    """



    """
            measured CCS vs. TJM CCS
    """
    """
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111)


    ms, mew = 7, 1.25
    ax.plot([100, 250], [100, 250], lin(n, *mopt), ls='--', c='grey', linewidth=1.5)
    ax.plot(bacs[1], tjm[1], 'ko', ms=ms, label='parents', mew=mew)
    ax.plot(bacs_O[1], tjm_O[1], 'b^', ms=ms, label='+O metab.', mew=mew)

    ax.set_xlim([140, 220])
    ax.set_ylim([140, 220])

    ax.set_xticks([140, 160, 180, 200, 220])
    ax.set_yticks([140, 160, 180, 200, 220])

    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.set_ylabel(r'theo. CCS ($\AA^2$)')
    ax.set_xlabel(r'exp. CCS ($\AA^2$)')
    ax.legend(fontsize=10)
    plt.tight_layout()

    # save figure
    png = 'BACs_ccsvsccs_comp.png'
    plt.savefig(png, dpi=400, bbox_inches='tight')

    #plt.show()
    plt.close()
    """


if __name__ == '__main__':
    main()



