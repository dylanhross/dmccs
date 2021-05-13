"""  """

from util import get_ccs_dist, plot_replicate_ccs_dists, two_peak_fit


# raw files
rf1 = 'D:/DMIM_HT/data/20191025/DHR20191025_p721_I10_pos.raw'
rf2 = 'D:/DMIM_HT/data/20191216/DHR20191216_p721_I10_pos.raw'
rf3 = 'D:/DMIM_HT/data/20191217/DHR20191217_p721_I10_pos.raw'

# ccs calibration files
cf1 = 'D:/DMIM_HT/analysis/20191025/combined_cal.pickle'
cf2 = 'D:/DMIM_HT/analysis/20191216/combined_cal.pickle'
cf3 = 'D:/DMIM_HT/analysis/20191217/combined_cal.pickle'

# get CCS distributions
cd1 = get_ccs_dist(rf1, 2, 352.1472, 0.05, cf1)
cd2 = get_ccs_dist(rf2, 1, 352.1472, 0.05, cf2)
cd3 = get_ccs_dist(rf3, 1, 352.1472, 0.05, cf3)

# plot atds 
plot_replicate_ccs_dists('lome', cd1, cd2, cd3)
plot_replicate_ccs_dists('lome', cd1, cd2, cd3, norm=True, x_range=[140, 220])

# fit CCS distributions to get CCS values
print()
ccs1, ccs2 = two_peak_fit(cd1[1], cd1[2], 165., 180., 1.)
print('CCS 1: {:.2f} Ang^2, CCS 2: {:.2f} Ang^2'.format(ccs1, ccs2))
ccs1, ccs2 = two_peak_fit(cd2[1], cd2[2], 165., 180., 1.)
print('CCS 1: {:.2f} Ang^2, CCS 2: {:.2f} Ang^2'.format(ccs1, ccs2))
ccs1, ccs2 = two_peak_fit(cd3[1], cd3[2], 165., 180., 1.)
print('CCS 1: {:.2f} Ang^2, CCS 2: {:.2f} Ang^2'.format(ccs1, ccs2))
print()

