#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3
"""
"""


from matplotlib import pyplot as plt
import json
import numpy as np
from scipy.optimize import curve_fit
import sys


fset = sys.argv[1]


# set up plot fonts
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica', 'Arial']
rcParams['font.size'] = 11


fig = plt.figure(figsize=(3, 4))
ax = fig.add_subplot(111)

# common settings for the same style across plots
bs = {
    'linewidth': 1.1, 'align': 'center', 'width': 0.8, 'capstyle': 'round', 'capsize': 5, 
    'error_kw': {
        'elinewidth': 1.1, 'ecolor': 'k'
    }
}


y = {
    'CUST': [2.210526517599249985e+02, 2.259224234838216319e+02],
    'MQN': [2.134326154759856990e+02, 2.134326154759856990e+02],
    'MD3D': [2.235046941620169036e+02, 2.179685099344615935e+02],
    'COMB': [1.938115028340808408e+02, 1.934436461179425635e+02]
}


c = 'k'
ax.bar([1], [229.5], edgecolor=c, color=c, fill=False, hatch='///', **bs)
ax.bar([1], [229.5], edgecolor=c, ecolor=c, fill=False, **bs)
c = 'b'
ax.bar([2], [235.9], edgecolor=c, color=c, fill=False, hatch='///', **bs)
ax.bar([2], [235.9], edgecolor=c, ecolor=c, fill=False, **bs)
c = 'k'
ax.bar([3], [y[fset][0]], edgecolor=c, color=c, **bs)
ax.bar([3], [y[fset][0]], edgecolor=c, ecolor=c, fill=False, **bs)
c = 'b'
ax.bar([4], [y[fset][1]], edgecolor=c, color=c, **bs)
ax.bar([4], [y[fset][1]], ecolor=c, edgecolor=c, fill=False, **bs)

ax.set_xticks([1, 2, 3, 4])
ax.set_xticklabels(['prot. A (exp.)', 'prot. B (exp.)', 'prot. A (pred.)', 'prot. B (pred.)'], fontstyle='italic', fontsize=10, rotation='vertical')


ax.set_ylim([190, 240])


for d in ['top', 'right']:
    ax.spines[d].set_visible(False)
ax.set_ylabel(r'CCS ($\AA^2$)')


plt.tight_layout()

# save figure
png = 'cefprox_{}_modeling_ccscomp.png'.format(fset)
plt.savefig(png, dpi=400, bbox_inches='tight')

#plt.show()
plt.close()





