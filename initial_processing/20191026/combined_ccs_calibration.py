"""
"""


from pandas import read_excel
import pickle

from dhrmasslynxapi.ccs_calibration import CCSCalibrationRaw



# get the calibrant info (m/z and lit. CCS) from the excel file
xl = read_excel('calibration.xlsx', sheet_name=None, index_col=None, header=0)
polya_mz, polya_ccs =  xl['polyala'].values.T.tolist()
drugs_mz, drugs_ccs =  xl['drugs'].values.T.tolist()

# calibration data file pairs
polya_base = 'F:/DMIM_HT/data/20191026/DHR20191026_polyala_pos_{:1d}.raw'
drugs_base = 'F:/DMIM_HT/data/20191026/DHR20191026_drugs_pos_{:1d}.raw'
cal_data_pairs = [[polya_base.format(i), drugs_base.format(i)] for i in range(1, 6)]

# set up a list of CCS calibrations based on each separate pair of calibration data files
individual_cals = [CCSCalibrationRaw(pair, [polya_mz, drugs_mz], [polya_ccs, drugs_ccs], mass_window=0.02) for pair in cal_data_pairs]

# dump the individual calibration curve figures
for i in range(5):
    individual_cals[i].cal_curve_figure('individual_cal_{:1d}.png'.format(i + 1))

# prepare a combined CCS calibration using all of the individual pairs
comb_dfiles, comb_mz, comb_ccs = [], [], []
for pair in cal_data_pairs:
    comb_dfiles.append(pair[0])
    comb_mz.append(polya_mz)
    comb_ccs.append(polya_ccs)
    comb_dfiles.append(pair[1])
    comb_mz.append(drugs_mz)
    comb_ccs.append(drugs_ccs)
combined_cal = CCSCalibrationRaw(comb_dfiles, comb_mz, comb_ccs, mass_window=0.05)
combined_cal.cal_curve_figure('combined_cal.png')

# save the calibration for later use
with open('combined_cal.pickle', 'wb') as pfile:
    pickle.dump(combined_cal, pfile)
