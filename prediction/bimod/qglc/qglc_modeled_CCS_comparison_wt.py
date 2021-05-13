#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3
"""
"""


from matplotlib import pyplot as plt
import numpy as np
import sys


fset = sys.argv[1]

# set up plot fonts
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica', 'Arial']
rcParams['font.size'] = 11


# common settings for the same style across plots
f_size = (3.5, 4)
bs = {
    'linewidth': 1., 'align': 'center', 'width': 0.75, 'capstyle': 'round', 'capsize': 6, 
    'error_kw': {
        'elinewidth': 1., 'ecolor': 'k'
    }
}
"""
bs = {
    'fill': False, 'linewidth': 2, 'align': 'center', 'width': 0.8, 'capstyle': 'round', 'capsize': 6, 'hatch': '//'
}
"""


fig = plt.figure(figsize=f_size)

ax = fig.add_subplot(111)
ax.axhline(202.8, c='k', ls=':', lw=1.5)
ax.axhline(209.5, c='k', ls=':', lw=1.5)

labels = ['7', '5', '3', "3'", "4'"]
x = [3, 4, 5, 2, 1]
ccs = np.loadtxt('y_{}.txt'.format(fset))
ec = ['orchid', 'yellow', 'lightpink', 'royalblue', 'darkorange']

for x_, ccs_, ec_ in zip(x, ccs, ec):
    ax.bar(x_, ccs_, edgecolor=ec_, ecolor=ec_, color=ec_, fill=True, **bs)
    ax.bar(x_, ccs_, fill=False, **bs)

ax.set_xticks([1, 2, 3, 4, 5])
ax.set_xticklabels(labels, fontstyle='italic', fontsize=14)



#ax.set_xlim([150, 400])
ax.set_ylim([190, 215])


for d in ['top', 'right']:
    ax.spines[d].set_visible(False)
ax.set_ylabel(r'CCS ($\AA^2$)')
#ax.set_xlabel('m/z')

# save figure
png = 'qglc_{}_modeled_comp_wt.png'.format(fset)
plt.savefig(png, dpi=400, bbox_inches='tight')

plt.tight_layout()
#plt.show()
#plt.close()

