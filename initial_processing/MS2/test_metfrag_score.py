""" """


import numpy as np

from metfrag import fragmenter_score


mz, i = np.loadtxt('cipr_MS2_peak2_nofilter.txt', unpack=True)
md1 = {
    'adduct': '[M+H]+',
    'formula': 'C17H18FN3O3',
    'monoiso': 331.1332,
    'inchi': 'InChI=1S/C17H18FN3O3/c18-13-7-11-14(8-15(13)20-5-3-19-4-6-20)21(10-1-2-10)9-12(16(11)22)17(23)24/h7-10,19H,1-6H2,(H,23,24)',
    'inchi_key': 'MYSWGUAQZAJSOK-UHFFFAOYSA-N'
}
md2 = {
    'adduct': '[M+H]+',
    'formula': 'C18H26ClN3O',
    'monoiso': 335.1764,
    'inchi': 'InChI=1S/C18H26ClN3O/c1-3-22(11-12-23)10-4-5-14(2)21-17-8-9-20-18-13-15(19)6-7-16(17)18/h6-9,13-14,23H,3-5,10-12H2,1-2H3,(H,20,21)',
    'inchi_key': 'MXXSMGPRMXLTPCZ-UHFFFAOYSA-N'
}
fscore1 = fragmenter_score(mz, i, 1000, md1)
fscore2 = fragmenter_score(mz, i, 1000, md2)
print(fscore1, fscore2)


# set up reader
#rdr = MassLynxReader(rawf)

# get the spectra
#p1_m, p1_i = rdr.get_spectrum(3, peak1_dt - dt_tol, peak1_dt + dt_tol)

