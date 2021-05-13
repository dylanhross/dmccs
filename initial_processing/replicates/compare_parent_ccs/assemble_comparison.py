"""

"""


from sqlite3 import connect
import numpy as np
from matplotlib import pyplot as plt

from dmim_analysis.util import remove_counter_ions


# set up database connections 
con_c3sdb = connect('C3S.db')
cur_c3sdb = con_c3sdb.cursor()
con_dmim = connect('replicate_data.db')
cur_dmim = con_dmim.cursor()


# set up the queries to gather the necessary data
qry_c3sdb = 'SELECT name, adduct, mz, ccs FROM master WHERE src_tag="hine0817";'
qry_dmim = """
SELECT cmpd, mz, ccs_avg, ccs_sd FROM plate_1 WHERE ccs_reps>1 AND cmpd NOT LIKE "%_met%" 
    UNION ALL SELECT cmpd, mz, ccs_avg, ccs_sd FROM plate_2 WHERE ccs_reps>1 AND cmpd NOT LIKE "%_met%"
    UNION ALL SELECT cmpd, mz, ccs_avg, ccs_sd FROM plate_3 WHERE ccs_reps>1 AND cmpd NOT LIKE "%_met%"
    UNION ALL SELECT cmpd, mz, ccs_avg, ccs_sd FROM plate_4 WHERE ccs_reps>1 AND cmpd NOT LIKE "%_met%"
    UNION ALL SELECT cmpd, mz, ccs_avg, ccs_sd FROM plate_5 WHERE ccs_reps>1 AND cmpd NOT LIKE "%_met%"
    UNION ALL SELECT cmpd, mz, ccs_avg, ccs_sd FROM plate_6 WHERE ccs_reps>1 AND cmpd NOT LIKE "%_met%"
    UNION ALL SELECT cmpd, mz, ccs_avg, ccs_sd FROM plate_7 WHERE ccs_reps>1 AND cmpd NOT LIKE "%_met%"
    UNION ALL SELECT cmpd, mz, ccs_r1, ccs_sd FROM plate_1 WHERE ccs_reps=1 AND cmpd NOT LIKE "%_met%"
    UNION ALL SELECT cmpd, mz, ccs_r1, ccs_sd FROM plate_2 WHERE ccs_reps=1 AND cmpd NOT LIKE "%_met%"
    UNION ALL SELECT cmpd, mz, ccs_r1, ccs_sd FROM plate_3 WHERE ccs_reps=1 AND cmpd NOT LIKE "%_met%"
    UNION ALL SELECT cmpd, mz, ccs_r1, ccs_sd FROM plate_4 WHERE ccs_reps=1 AND cmpd NOT LIKE "%_met%"
    UNION ALL SELECT cmpd, mz, ccs_r1, ccs_sd FROM plate_5 WHERE ccs_reps=1 AND cmpd NOT LIKE "%_met%"
    UNION ALL SELECT cmpd, mz, ccs_r1, ccs_sd FROM plate_6 WHERE ccs_reps=1 AND cmpd NOT LIKE "%_met%"
    UNION ALL SELECT cmpd, mz, ccs_r1, ccs_sd FROM plate_7 WHERE ccs_reps=1 AND cmpd NOT LIKE "%_met%"
"""


# fetch data from the databases
cmpds_c3sdb, mzs_c3sdb, ccss_c3sdb = [], [], []
for name, adduct, mz, ccs in cur_c3sdb.execute(qry_c3sdb).fetchall():
    # concatenate the compound name and adduct into one identifier
    cmpds_c3sdb.append('{}_{}'.format(remove_counter_ions(name), adduct))
    mzs_c3sdb.append(float(mz))
    ccss_c3sdb.append(float(ccs))

cmpds_dmim, mzs_dmim, ccss_dmim, ccssds_dmim = [], [], [], []
for cmpd, mz, ccs, ccs_sd in cur_dmim.execute(qry_dmim).fetchall():
    cmpds_dmim.append(cmpd)
    mzs_dmim.append(float(mz))
    ccss_dmim.append(float(ccs))
    ccssds_dmim.append(float(ccs_sd) if ccs_sd is not None else None)


def percent_error(y_true, y_new):
    return 100. * (y_new - y_true) / y_true


def find_matching(compound):
    for cmpd, mz, ccs, ccssd in zip(cmpds_dmim, mzs_dmim, ccss_dmim, ccssds_dmim):
        if compound == cmpd or compound.replace(',', '-') == cmpd:
            return mz, ccs, ccssd
    return None, None, None


count = 0
ccs_pe_cmpds = []
ccs_pe = []
for cmpd, mz, ccs in zip(cmpds_c3sdb, mzs_c3sdb, ccss_c3sdb):
    print('{:50s}'.format(cmpd), end=' ... ')
    mz2, ccs2, ccssd = find_matching(cmpd)
    if mz2 is not None:
        count += 1
        s = 'MZ {:9.4f} {:9.4f} ({:+7d} ppm) CCS {:6.2f} {:6.2f} ({:+6.2f} %)'
        print(s.format(mz, mz2, int(10000. * percent_error(mz, mz2)), ccs, ccs2, percent_error(ccs, ccs2)))
    
        ppm = 10000. * percent_error(mz, mz2)
        if np.abs(ppm) < 10000:
            ccs_pe.append(percent_error(ccs, ccs2))
            ccs_pe_cmpds.append(cmpd)
    else:
        print('NOT FOUND')


print('\n\n\n')
print(count)
print(len(ccs_pe))
print('mean {:.2f} %'.format(np.mean(ccs_pe)))
print('median {:.2f} %'.format(np.median(ccs_pe)))
print('<=5%: {}'.format(len(np.array(ccs_pe)[np.abs(ccs_pe) <= 5])))
print('<=3%: {}'.format(len(np.array(ccs_pe)[np.abs(ccs_pe) <= 3])))
print('<=1%: {}'.format(len(np.array(ccs_pe)[np.abs(ccs_pe) <= 1])))


bins = np.linspace(-6, 6, 61).tolist()
plt.hist(ccs_pe, bins=bins, color='b', edgecolor='k')
plt.xlabel('% error')
plt.ylabel('count')
plt.show()


