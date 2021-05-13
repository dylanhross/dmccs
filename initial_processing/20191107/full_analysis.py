"""
"""


import os
from csv import reader
from pickle import load, dump

from dmim_analysis.analysis import analyze_compound
from dmim_analysis.plots import parent_atd_with_fit_and_structure, metabolite_atd_with_fit_and_structure


with open('combined_cal.pickle', 'rb') as pf:
    cal = load(pf)


i = 1
with open('plate_4.csv', 'r') as f:
    rdr = reader(f)
    for name, f_rxn, f_ctrl in rdr:
        print('analyzing {} ...'.format(name), end=' ')
        f_rxn = 'D:/DMIM_HT/data/20191107/' + f_rxn
        f_ctrl = 'D:/DMIM_HT/data/20191107/' + f_ctrl
        try:
            cmpd = analyze_compound(i, name, f_rxn, 0, 2, 0.05, 'pos', cal, check_ctrl=True, ctrl_raw_file=f_ctrl)
        except KeyboardInterrupt:
            print('\n! EXITING !')
            exit()
        except Exception as e:
            print(e)
            cmpd = None
        if cmpd:
            save_name = name.replace(' ', '_').replace(':', '.')
            prefix = 'compounds/{}'.format(save_name)
            if not os.path.isdir(prefix):
                os.mkdir(prefix)
            with open(os.path.join(prefix, '{}.pickle'.format(save_name, save_name)), 'wb') as pf:
                dump(cmpd, pf)
            try:
                parent_atd_with_fit_and_structure(cmpd, prefix)
                metabolite_atd_with_fit_and_structure(cmpd, prefix)
            except Exception as e:
                print(e)
                print("unable to generate plot.")
            print('ok')
        else:
            print('failed')
        i += 1



