"""  """

from util import (
    get_ccs_dist, plot_replicate_ccs_dists, plot_replicate_atds, two_peak_fit, plot_fragment_spectra,
    plot_fragment_atd
)


# raw files
rf1 = 'D:/DMIM_HT/data/20200127/DHR20200127_p723_O07_pos.raw'
rf2 = 'D:/DMIM_HT/data/20191209/DHR20191209_p723_O07_pos.raw'
rf3 = 'D:/DMIM_HT/data/20191210/DHR20191210_p723_O07_pos.raw'

# ccs calibration files
cf1 = 'D:/DMIM_HT/analysis/20200127/combined_cal.pickle'
cf2 = 'D:/DMIM_HT/analysis/20191209/combined_cal.pickle'
cf3 = 'D:/DMIM_HT/analysis/20191210/combined_cal.pickle'

# get CCS distributions
cd1 = get_ccs_dist(rf1, 2, 321.1363, 0.05, cf1)
cd2 = get_ccs_dist(rf2, 1, 321.1363, 0.05, cf2)
cd3 = get_ccs_dist(rf3, 1, 321.1363, 0.05, cf3)

# plot atds 
plot_replicate_ccs_dists('enox', cd1, cd2, cd3)
plot_replicate_ccs_dists('enox', cd1, cd2, cd3, norm=True, x_range=[140, 220])
plot_replicate_atds('enox', cd1, cd2, cd3)
plot_replicate_atds('enox', cd1, cd2, cd3, norm=True)

# fit ATDs to get drift times
print()
dt1, dt2 = two_peak_fit(cd1[0], cd1[2], 3.3, 3.8, 0.1)
print('dt 1: {:.2f} ms, dt 2: {:.2f} ms'.format(dt1, dt2))
dt1, dt2 = two_peak_fit(cd2[0], cd2[2], 3.3, 3.8, 0.1)
print('dt 1: {:.2f} ms, dt 2: {:.2f} ms'.format(dt1, dt2))
dt1, dt2 = two_peak_fit(cd3[0], cd3[2], 3.3, 3.8, 0.1)
print('dt 1: {:.2f} ms, dt 2: {:.2f} ms'.format(dt1, dt2))

# fit CCS distributions to get CCS values
print()
ccs1, ccs2 = two_peak_fit(cd1[1], cd1[2], 165., 180., 1.)
print('CCS 1: {:.2f} Ang^2, CCS 2: {:.2f} Ang^2'.format(ccs1, ccs2))
ccs1, ccs2 = two_peak_fit(cd2[1], cd2[2], 165., 180., 1.)
print('CCS 1: {:.2f} Ang^2, CCS 2: {:.2f} Ang^2'.format(ccs1, ccs2))
ccs1, ccs2 = two_peak_fit(cd3[1], cd3[2], 165., 180., 1.)
print('CCS 1: {:.2f} Ang^2, CCS 2: {:.2f} Ang^2'.format(ccs1, ccs2))
print()

# MS/MS spectra for each peak
plot_fragment_spectra('enox', rf1, 3.2, 3.7, 0.1, 50, 350)
plot_fragment_spectra('enox_fA', rf1, 3.2, 3.7, 0.1, 292.8, 293.6)
plot_fragment_spectra('enox_fB', rf1, 3.2, 3.7, 0.1, 257.8, 258.6)
plot_fragment_spectra('enox_fC', rf1, 3.2, 3.7, 0.1, 233.8, 234.6)
plot_fragment_spectra('enox_fD', rf1, 3.2, 3.7, 0.1, 231.8, 232.6)
plot_fragment_spectra('enox_fE', rf1, 3.2, 3.7, 0.1, 229.8, 230.6)
plot_fragment_spectra('enox_fF', rf1, 3.2, 3.7, 0.1, 203.8, 204.6)
plot_fragment_spectra('enox_fG', rf1, 3.2, 3.7, 0.1, 83.8, 84.6)

# plot ATDs for each fragment
plot_fragment_atd('enox_fA', rf1, 3, 293.15, 0.05, cd1, dt_scale=1.08)
plot_fragment_atd('enox_fB', rf1, 3, 258.15, 0.05, cd1, dt_scale=1.08)
plot_fragment_atd('enox_fC', rf1, 3, 234.12, 0.05, cd1, dt_scale=1.08)
plot_fragment_atd('enox_fD', rf1, 3, 232.09, 0.05, cd1, dt_scale=1.08)
plot_fragment_atd('enox_fE', rf1, 3, 230.09, 0.05, cd1, dt_scale=1.08)
plot_fragment_atd('enox_fF', rf1, 3, 204.06, 0.05, cd1, dt_scale=1.08)
plot_fragment_atd('enox_fG', rf1, 3, 84.09, 0.05, cd1, dt_scale=1.08)

