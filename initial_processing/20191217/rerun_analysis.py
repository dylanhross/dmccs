"""
"""


import os
from csv import reader
from pickle import load, dump

from dmim_analysis.analysis import analyze_compound_targeted
from dmim_analysis.plots import atd_with_fit


with open('combined_cal.pickle', 'rb') as pf:
    cal = load(pf)


with open('D:/DMIM_HT/analysis/replicates/rerun_p1r3.csv', 'r') as f:
    rdr = reader(f)
    for name_adduct, well_loc, mz, ccs_rep1 in rdr:
        mz, ccs_rep1 = float(mz), float(ccs_rep1)
        print('analyzing {} ...'.format(name_adduct), end=' ')
        f_rxn = 'D:/DMIM_HT/data/20191217/DHR20191217_p721_{}_pos.raw'.format(well_loc)
        try:
            cmpd = analyze_compound_targeted(name_adduct, mz, f_rxn, 0, 1, 0.0825, cal)
        except KeyboardInterrupt:
            print('\n! EXITING !')
            exit()
        except Exception as e:
            print(e)
            cmpd = None
        if cmpd:
            p_err = 100. * (cmpd['ccs'] - ccs_rep1) / ccs_rep1
            print('CCS: {:.2f} <-> rep1 CCS: {:.2f} ({:+.1f} %) ...'.format(cmpd['ccs'], ccs_rep1, p_err), end=' ')
            save_name = name_adduct.replace(' ', '_')
            prefix = 'compounds/'
            with open(os.path.join(prefix, '{}.pickle'.format(save_name)), 'wb') as pf:
                dump(cmpd, pf)
            atd_with_fit(name_adduct, mz, cmpd['ccs'], cmpd['meta']['atd_fit_params'], cmpd['meta']['atd'], 
                         cmpd['meta']['lc_fit_params'], cmpd['meta']['lc'], 0.0825, prefix)
            print('ok')
        else:
            print('failed')




