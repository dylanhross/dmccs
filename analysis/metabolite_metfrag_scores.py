#!/usr/bin/python3

from matplotlib import pyplot as plt
from numpy import loadtxt, array, log10
from numpy.random import normal
from sqlite3 import connect

from matplotlib import rcParams
rcParams['font.size'] = 8


# fetch scores
qry_metabolite_metfrag_scores = """
SELECT 
    metfrag_score
FROM
    plate_{n}_ms2
JOIN 
    plate_{n}
    ON plate_{n}.dmim_id = plate_{n}_ms2.dmim_id
WHERE
    -- only metabolites
    met_n > 0
    AND metfrag_score > 0
;"""
con = connect('DMIM_v1.1.db')
cur = con.cursor()
scores = []
for n in range(1, 8):
    for score in cur.execute(qry_metabolite_metfrag_scores.format(n=n)).fetchall():
        scores.append(score)
metfrag_scores = array(scores)
con.close()
#metfrag_scores = loadtxt('build/DMIM_metfrag_scores_prefilter.txt')

print(len(metfrag_scores))


# plot distributions of the metfrag scores for top-10 ranked compounds
xs = normal(loc=1, scale=0.1, size=len(metfrag_scores))
fig = plt.figure(figsize=(2, 3))
ax = fig.add_subplot(111)
ax.scatter(xs, log10(metfrag_scores), marker='.', s=8, c='r', edgecolors='none', alpha=0.4)
flierprops = dict(marker='o', markeredgecolor='k', markersize=2, linestyle='none', markeredgewidth=0.5, fillstyle='none')
bpd = ax.boxplot(log10(metfrag_scores), flierprops=flierprops)
for l in bpd['medians']:
    l.set_color('k')
ax.set_ylabel(r'$log_{10}$(MetFrag score)')
for d in ['top', 'right']:
    ax.spines[d].set_visible(False)
ax.set_xticks([])
ax.set_xlabel('metabolites')
ax.set_xlim([0.5, 1.5])
plt.savefig('DMIM_metabolites_mf_scores.png', dpi=350, bbox_inches='tight')
plt.close()
