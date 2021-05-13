#!/usr/local/Cellar/python@3.9/3.9.1_6/bin/python3 
""" """

from matplotlib import pyplot as plt, rcParams
from numpy import mean, std, loadtxt
from numpy.random import normal
import sys


# set the global fontsize on plots
rcParams['font.size'] = 8

# hard code the measured CCS values
meas_ccs = {
    'cipr_B': [172.67, 172.07, 171.69],
    'cipr_A': [183.93, 184.11, 183.94],
    'enox_B': [167.74, 167.29, 167.34],
    'enox_A': [183.43, 183.35, 183.44],
    'enro_B': [158.30, 157.65],
    'enro_A': [188.32, 188.32],
    'levo_B': [153.63, 178.28, 178.21],
    'levo_A': [183.80, 187.72, 187.45],
    'lome_B': [175.16, 174.73, 174.76],
    'lome_A': [187.85, 187.06, 187.16],
    'norf_B': [169.43, 169.43, 169.60],
    'norf_A': [182.93, 182.63, 182.84],
    'orbi_B': [184.55, 184.28],
    'orbi_A': [196.13, 195.78],
    'pefl_B': [164.75, 165.28, 165.52],
    'pefl_A': [179.59, 174.60, 174.17]
}


def fetch_data(label, fset, n):
    """ """
    meas_A = meas_ccs['{}_A'.format(label)]
    meas_B = meas_ccs['{}_B'.format(label)]
    y = loadtxt('{}_y_{}.txt'.format(label, fset))
    calc_A = y[:n]
    calc_B = y[n:2 * n]
    calc_C = y[2 * n:3 * n]
    calc_D = y[3 * n:]
    return meas_A, meas_B, calc_A, calc_B, calc_C, calc_D


def single_plot(label, fset, meas_A, meas_B, calc_A, calc_B, calc_C, calc_D, ymin):
    """ """
    fig = plt.figure(figsize=(1.5, 2.))
    ax = fig.add_subplot(111)
    
    width = 0.4
    x = [1 - width / 2, 1 + width / 2, 2 - width / 2, 2 + width / 2, 2 + (3 / 2) * width, 2 + (5 / 2) * width]
    y = [meas_A, meas_B, calc_A, calc_B, calc_C, calc_D]
    height = [mean(_) - ymin for _ in y]
    yerr = [std(meas_A), std(meas_B), std(calc_A), std(calc_B), std(calc_C), std(calc_D)]
    color = ['#333333', '#999999', '#0A2F51', '#137177', '#1D9A6C', '#74C67A']
    ax.bar(x, height, width=width, yerr=yerr, color=color, bottom=ymin)
    # add scatter dots for all the points
    xx = []
    yy = []
    for i in range(6):
        for j in range(len(y[i])):
            xx.append(x[i] + normal(loc=0, scale=0.05))
            yy.append(y[i][j])
    ax.plot(xx, yy, 'ko', ms=2, fillstyle='none', mew=0.5)
    
    #ax.set_xticks(x)
    #ax.set_xticklabels(['meas. A', 'meas. B', 'pred. A', 'pred. B'], rotation=45)
    ax.set_xticks([1, 2. + width])
    ax.set_xticklabels(['meas.', 'pred.'])

    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)

    ax.set_ylabel(r'CCS ($\AA^2$)', fontsize=8)
    fig.savefig('{}_{}_meas_pred_ccs.png'.format(label, fset), dpi=400, bbox_inches='tight')
    plt.close()


def main():
    """ main execution sequence """
    n = int(sys.argv[1])
    label = sys.argv[2]
    ymin = float(sys.argv[3])

    meas_A, meas_B, calc_A, calc_B, calc_C, calc_D = fetch_data(label, 'CUST', n)
    single_plot(label, 'CUST', meas_A, meas_B, calc_A, calc_B, calc_C, calc_D, ymin)
    meas_A, meas_B, calc_A, calc_B, calc_C, calc_D = fetch_data(label, 'MQN', n)
    single_plot(label, 'MQN', meas_A, meas_B, calc_A, calc_B, calc_C, calc_D, ymin)
    meas_A, meas_B, calc_A, calc_B, calc_C, calc_D = fetch_data(label, 'MD3D', n)
    single_plot(label, 'MD3D', meas_A, meas_B, calc_A, calc_B, calc_C, calc_D, ymin)
    meas_A, meas_B, calc_A, calc_B, calc_C, calc_D = fetch_data(label, 'COMB', n)
    single_plot(label, 'COMB', meas_A, meas_B, calc_A, calc_B, calc_C, calc_D, ymin)


if __name__ == '__main__':
    main()


