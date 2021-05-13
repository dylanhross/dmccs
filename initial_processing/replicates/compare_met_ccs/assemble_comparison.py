"""

"""


from sqlite3 import connect
import numpy as np
from matplotlib import pyplot as plt

from dmim_analysis.util import remove_counter_ions


# set up database connections 
con_dmim = connect('replicate_data.db')
cur_dmim = con_dmim.cursor()


# set up the queries to gather the necessary data
qry_dmim = """
SELECT cmpd, mz, ccs_avg, ccs_sd FROM plate_1 WHERE ccs_reps>1 AND cmpd LIKE "%_met%" 
    UNION ALL SELECT cmpd, mz, ccs_avg, ccs_sd FROM plate_2 WHERE ccs_reps>1 AND cmpd LIKE "%_met%"
    UNION ALL SELECT cmpd, mz, ccs_avg, ccs_sd FROM plate_3 WHERE ccs_reps>1 AND cmpd LIKE "%_met%"
    UNION ALL SELECT cmpd, mz, ccs_avg, ccs_sd FROM plate_4 WHERE ccs_reps>1 AND cmpd LIKE "%_met%"
    UNION ALL SELECT cmpd, mz, ccs_avg, ccs_sd FROM plate_5 WHERE ccs_reps>1 AND cmpd LIKE "%_met%"
    UNION ALL SELECT cmpd, mz, ccs_avg, ccs_sd FROM plate_6 WHERE ccs_reps>1 AND cmpd LIKE "%_met%"
    UNION ALL SELECT cmpd, mz, ccs_avg, ccs_sd FROM plate_7 WHERE ccs_reps>1 AND cmpd LIKE "%_met%"
"""


cmpds_dmim, mzs_dmim, ccss_dmim, ccssds_dmim, ccs_rsd = [], [], [], [], []
for cmpd, mz, ccs, ccs_sd in cur_dmim.execute(qry_dmim).fetchall():
    cmpds_dmim.append(cmpd)
    mzs_dmim.append(float(mz))
    ccss_dmim.append(float(ccs))
    ccssds_dmim.append(float(ccs_sd))
    ccs_rsd.append((100. * ccs_sd) / ccs)

for cmpd, mz, ccs, rsd in zip(cmpds_dmim, mzs_dmim, ccss_dmim, ccs_rsd):
    print('{:60s}'.format(cmpd), end=' ... ')
    s = 'MZ {:9.4f} CCS {:6.2f} ({:5.2f} %)'
    print(s.format(mz, ccs, rsd))


print('\n\n\n')
print(len(ccs_rsd))
print('mean {:.2f} %'.format(np.mean(ccs_rsd)))
print('median {:.2f} %'.format(np.median(ccs_rsd)))


bins = np.linspace(0, 2, 21).tolist()
plt.hist(ccs_rsd, bins=bins, color='b', edgecolor='k')
plt.xlabel('RSD %')
plt.ylabel('count')
plt.show()


