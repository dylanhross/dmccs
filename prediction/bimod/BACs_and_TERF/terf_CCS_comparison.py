#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3
"""
"""


from matplotlib import pyplot as plt
import json
import numpy as np
from scipy.optimize import curve_fit
import sys

from BACs_ccs_vs_mz import compaction_factor


fset = sys.argv[1]

# set up plot fonts
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica', 'Arial']
rcParams['font.size'] = 11

# compute weighted average of predicted CCS
a = np.loadtxt('TERF_y_{}.txt'.format(fset))
b = np.loadtxt('TERF_cluster_sizewt.txt', usecols=(1))
c = np.loadtxt('TERFO_y_{}.txt'.format(fset))
d = np.loadtxt('TERFO_cluster_sizewt.txt', usecols=(1))
wavg = sum(a * b) / sum(b)
wavgO = sum(c * d) / sum(d)

fig = plt.figure(figsize=(3, 4))
ax = fig.add_subplot(111)

# common settings for the same style across plots
bs = {
    'linewidth': 1.1, 'align': 'center', 'width': 0.8, 'capstyle': 'round', 'capsize': 5, 
    'error_kw': {
        'elinewidth': 1.1, 'ecolor': 'k'
    }
}

c = 'k'
ax.bar([1], [227.04], edgecolor=c, color=c, fill=False, hatch='///', **bs)
ax.bar([1], [227.04], edgecolor=c, ecolor=c, fill=False, **bs)
c = 'b'
ax.bar([2], [222.82], edgecolor=c, color=c, fill=False, hatch='///', **bs)
ax.bar([2], [222.82], edgecolor=c, ecolor=c, fill=False, **bs)
ax.text(1.7, 225, '{:.3f}'.format(compaction_factor(472.3376, 227.04, 488.3268, 222.82)), fontsize=10, c='b')
c = 'k'
ax.bar([3], [wavg], edgecolor=c, color=c, **bs)
ax.bar([3], [wavg], edgecolor=c, ecolor=c, fill=False, **bs)
c = 'b'
ax.bar([4], [wavgO], edgecolor=c, color=c, **bs)
ax.bar([4], [wavgO], ecolor=c, edgecolor=c, fill=False, **bs)
ax.text(3.7, wavgO + 2, '{:.3f}'.format(compaction_factor(472.3376, wavg, 488.3268, wavgO)), fontsize=10, c='b')

ax.set_xticks([1, 2, 3, 4])
ax.set_xticklabels(['TERF (exp.)', 'TERF+O (exp.)', 'TERF (pred.)', 'TERF+O (pred.)'], fontstyle='italic', fontsize=10, rotation='vertical')


ax.set_ylim([200, 245])


for d in ['top', 'right']:
    ax.spines[d].set_visible(False)
ax.set_ylabel(r'CCS ($\AA^2$)')


plt.tight_layout()

# save figure
png = 'terf_{}_modeling_ccscomp.png'.format(fset)
plt.savefig(png, dpi=400, bbox_inches='tight')

#plt.show()
plt.close()





