"""  """

from util import (
    get_ccs_dist, plot_replicate_ccs_dists, plot_replicate_atds, two_peak_fit, plot_fragment_spectra,
    plot_fragment_atd
)


# raw files
rf1 = 'D:/DMIM_HT/data/20191025/DHR20191025_p721_G17_pos.raw'
rf2 = 'D:/DMIM_HT/data/20191216/DHR20191216_p721_G17_pos.raw'
rf3 = 'D:/DMIM_HT/data/20191217/DHR20191217_p721_G17_pos.raw'

# ccs calibration files
cf1 = 'D:/DMIM_HT/analysis/20191025/combined_cal.pickle'
cf2 = 'D:/DMIM_HT/analysis/20191216/combined_cal.pickle'
cf3 = 'D:/DMIM_HT/analysis/20191217/combined_cal.pickle'

# get CCS distributions
cd1 = get_ccs_dist(rf1, 2, 320.1410, 0.05, cf1)
cd2 = get_ccs_dist(rf2, 1, 320.1410, 0.05, cf2)
cd3 = get_ccs_dist(rf3, 1, 320.1410, 0.05, cf3)

# plot atds 
plot_replicate_ccs_dists('norf', cd1, cd2, cd3)
plot_replicate_ccs_dists('norf', cd1, cd2, cd3, norm=True, x_range=[140, 220])
plot_replicate_atds('norf', cd1, cd2, cd3)
plot_replicate_atds('norf', cd1, cd2, cd3, norm=True)

# fit ATDs to get drift times
print()
dt1, dt2 = two_peak_fit(cd1[0], cd1[2], 2.8, 3.5, 0.1)
print('dt 1: {:.2f} ms, dt 2: {:.2f} ms'.format(dt1, dt2))
dt1, dt2 = two_peak_fit(cd2[0], cd2[2], 3.2, 3.7, 0.1)
print('dt 1: {:.2f} ms, dt 2: {:.2f} ms'.format(dt1, dt2))
dt1, dt2 = two_peak_fit(cd3[0], cd3[2], 3.2, 3.7, 0.1)
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
plot_fragment_spectra('norf', rf1, 2.9, 3.34, 0.1, 50, 350)
plot_fragment_spectra('norf_fA', rf1, 2.9, 3.34, 0.1, 301.8, 302.6)
plot_fragment_spectra('norf_fB', rf1, 2.9, 3.34, 0.1, 275.8, 276.6)
plot_fragment_spectra('norf_fC', rf1, 2.9, 3.34, 0.1, 255.8, 256.6)
plot_fragment_spectra('norf_fD', rf1, 2.9, 3.34, 0.1, 232.8, 233.6)

# plot ATDs for each fragment
plot_fragment_atd('norf_fA', rf1, 3, 301.13, 0.05, cd1, dt_scale=1.06)
plot_fragment_atd('norf_fB', rf1, 3, 276.16, 0.05, cd1, dt_scale=1.06)
plot_fragment_atd('norf_fC', rf1, 3, 256.16, 0.05, cd1, dt_scale=1.06)
plot_fragment_atd('norf_fD', rf1, 3, 233.11, 0.05, cd1, dt_scale=1.06)


