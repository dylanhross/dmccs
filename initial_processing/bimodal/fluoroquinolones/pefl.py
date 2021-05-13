"""  """

from util import get_ccs_dist, plot_replicate_ccs_dists, two_peak_fit


# raw files
rf1 = 'D:/DMIM_HT/data/20200127/DHR20200127_p723_L21_pos.raw'
rf2 = 'D:/DMIM_HT/data/20191209/DHR20191209_p723_L21_pos.raw'
rf3 = 'D:/DMIM_HT/data/20191210/DHR20191210_p723_L21_pos.raw'

# ccs calibration files
cf1 = 'D:/DMIM_HT/analysis/20200127/combined_cal.pickle'
cf2 = 'D:/DMIM_HT/analysis/20191209/combined_cal.pickle'
cf3 = 'D:/DMIM_HT/analysis/20191210/combined_cal.pickle'

# get CCS distributions
cd1 = get_ccs_dist(rf1, 2, 334.1567, 0.05, cf1)
cd2 = get_ccs_dist(rf2, 1, 334.1567, 0.05, cf2)
cd3 = get_ccs_dist(rf3, 1, 334.1567, 0.05, cf3)

# plot atds 
plot_replicate_ccs_dists('pefl', cd1, cd2, cd3)
plot_replicate_ccs_dists('pefl', cd1, cd2, cd3, norm=True, x_range=[140, 220])

# fit CCS distributions to get CCS values
print()
ccs1, ccs2 = two_peak_fit(cd1[1], cd1[2], 165., 180., 1.)
print('CCS 1: {:.2f} Ang^2, CCS 2: {:.2f} Ang^2'.format(ccs1, ccs2))
ccs1, ccs2 = two_peak_fit(cd2[1], cd2[2], 165., 180., 1.)
print('CCS 1: {:.2f} Ang^2, CCS 2: {:.2f} Ang^2'.format(ccs1, ccs2))
ccs1, ccs2 = two_peak_fit(cd3[1], cd3[2], 165., 180., 1.)
print('CCS 1: {:.2f} Ang^2, CCS 2: {:.2f} Ang^2'.format(ccs1, ccs2))
print()

