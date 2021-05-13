"""  """

from util import (
    get_ccs_dist, plot_replicate_ccs_dists, plot_replicate_atds, two_peak_fit, plot_fragment_spectra,
    plot_fragment_atd
)


# raw files
rf1 = 'D:/DMIM_HT/data/20200127/DHR20200127_p723_I16_pos.raw'
rf2 = 'D:/DMIM_HT/data/20191209/DHR20191209_p723_I16_pos.raw'
rf3 = 'D:/DMIM_HT/data/20191210/DHR20191210_p723_I16_pos.raw'

# ccs calibration files
cf1 = 'D:/DMIM_HT/analysis/20200127/combined_cal.pickle'
cf2 = 'D:/DMIM_HT/analysis/20191209/combined_cal.pickle'
cf3 = 'D:/DMIM_HT/analysis/20191210/combined_cal.pickle'

# get CCS distributions
cd1 = get_ccs_dist(rf1, 2, 334.1410, 0.05, cf1)
cd2 = get_ccs_dist(rf2, 1, 334.1410, 0.05, cf2)
cd3 = get_ccs_dist(rf3, 1, 334.1410, 0.05, cf3)

# plot atds and CCS distributions
plot_replicate_ccs_dists('cipr', cd1, cd2, cd3)
plot_replicate_ccs_dists('cipr', cd1, cd2, cd3, norm=True, x_range=[140, 220])
plot_replicate_atds('cipr', cd1, cd2, cd3)
plot_replicate_atds('cipr', cd1, cd2, cd3, norm=True)

# fit ATDs to get drift times
print()
dt1, dt2 = two_peak_fit(cd1[0], cd1[2], 3.4, 3.8, 0.1)
print('dt 1: {:.2f} ms, dt 2: {:.2f} ms'.format(dt1, dt2))
dt1, dt2 = two_peak_fit(cd2[0], cd2[2], 3.4, 3.8, 0.1)
print('dt 1: {:.2f} ms, dt 2: {:.2f} ms'.format(dt1, dt2))
dt1, dt2 = two_peak_fit(cd3[0], cd3[2], 3.4, 3.8, 0.1)
print('dt 1: {:.2f} ms, dt 2: {:.2f} ms'.format(dt1, dt2))

# fit CCS distributions to get CCS values
print()
ccs1, ccs2 = two_peak_fit(cd1[1], cd1[2], 170., 185., 1.)
print('CCS 1: {:.2f} Ang^2, CCS 2: {:.2f} Ang^2'.format(ccs1, ccs2))
ccs1, ccs2 = two_peak_fit(cd2[1], cd2[2], 170., 185., 1.)
print('CCS 1: {:.2f} Ang^2, CCS 2: {:.2f} Ang^2'.format(ccs1, ccs2))
ccs1, ccs2 = two_peak_fit(cd3[1], cd3[2], 170., 185., 1.)
print('CCS 1: {:.2f} Ang^2, CCS 2: {:.2f} Ang^2'.format(ccs1, ccs2))
print()

# MS/MS spectra for each peak
plot_fragment_spectra('cipr', rf1, 3.36, 3.74, 0.1, 50, 340, dump_txt=True)
plot_fragment_spectra('cipr_fA', rf1, 3.36, 3.74, 0.1, 313.8, 314.6)
plot_fragment_spectra('cipr_fB', rf1, 3.36, 3.74, 0.1, 287.8, 288.6)
plot_fragment_spectra('cipr_fC', rf1, 3.36, 3.74, 0.1, 270.8, 271.6)
plot_fragment_spectra('cipr_fD', rf1, 3.36, 3.74, 0.1, 244.8, 245.6)
plot_fragment_spectra('cipr_fE', rf1, 3.36, 3.74, 0.1, 230.8, 231.6)
plot_fragment_spectra('cipr_fF', rf1, 3.36, 3.74, 0.1, 202.8, 203.6)

# plot ATDs for each fragment
plot_fragment_atd('cipr_fA', rf1, 3, 314.15, 0.05, cd1, dt_scale=1.08)
plot_fragment_atd('cipr_fB', rf1, 3, 288.16, 0.05, cd1, dt_scale=1.08)
plot_fragment_atd('cipr_fC', rf1, 3, 271.03, 0.05, cd1, dt_scale=1.08)
plot_fragment_atd('cipr_fD', rf1, 3, 245.12, 0.05, cd1, dt_scale=1.08)
plot_fragment_atd('cipr_fE', rf1, 3, 231.06, 0.05, cd1, dt_scale=1.08)
plot_fragment_atd('cipr_fF', rf1, 3, 203.05, 0.05, cd1, dt_scale=1.08)
