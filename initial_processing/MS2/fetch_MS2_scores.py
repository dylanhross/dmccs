""" """

from sqlite3 import connect
import re
import pickle
import os
from glob import glob

from dmim_analysis.util import remove_counter_ions
from dhrmasslynxapi.reader import MassLynxReader


# set up db connection and cursors for running queries
con = connect('replicate_data_assigned_MS2.db')
cur1, cur2 = con.cursor(), con.cursor()


# add the MS2 tables to the database
for i in range(1, 8):
    cur1.execute("""CREATE TABLE plate_{n}_ms2 (cmpd TEXT NOT NULL, well TEXT NOT NULL, spectrum TEXT NOT NULL, adduct TEXT, monoiso REAL, formula TEXT, inchi TEXT, inchi_key TEXT, smi TEXT, metfrag_score REAL);""".format(n=i))


# define queries
qry1 = """SELECT cmpd, well FROM plate_{n}"""
qry2 = """INSERT INTO plate_{n}_ms2 (cmpd, well, spectrum, adduct, monoiso, smi) VALUES (?,?,?,?,?,?)"""
qry3 = """INSERT INTO plate_{n}_ms2 (cmpd, well, spectrum, adduct, monoiso, inchi) VALUES (?,?,?,?,?,?)"""


# define regex patterns
pat1 = re.compile(r'(.+)_(\[M(\+[HNaK]{1,2}(-H2O){0,1}){0,1}\]\+)')
pat2 = re.compile(r'(.*)_met([0-9]+)')


# returns the path to the pickle file with the compound data structure
def find_cmpd_pickle(name, plate_date, sub_commas=True):
    name = name.replace(' ', '_')
    for ppath in glob("D:/DMIM_HT/analysis/{}/compounds/*".format(plate_date)):
        dirname = os.path.split(ppath)[1]
        dirname_edit = remove_counter_ions(dirname, separator='_')
        if sub_commas:
            dirname_edit = dirname_edit.replace(',', '-')
        if name == dirname_edit:
            return ppath + '/' + dirname + '.pickle'
    return None


# returns a string containing the mass spectrum in plain-text format
def spectrum_as_txt(m, i, m_max):
    s = ''
    f = '{:.4f} {:d}\n'
    for m_, i_ in zip(m, i):
        if m_ > m_max:
            break
        s += f.format(m_, int(i_))
    return s


# iterate through all 7 plates
plate_dates = ['20191025', '20191026', '20200127', '20191107', '20191108', '20191113', '20191114']
for plate_n in range(1, 8):
    plate_date = plate_dates[plate_n - 1]
    # iterate through all compounds in the plate
    for cmpd, well in cur1.execute(qry1.format(n=plate_n)).fetchall():
        # parse compound name and adduct from the compound id
        mat1 = pat1.match(cmpd)
        name, adduct = mat1.group(1), mat1.group(2)
        # determine if this is a metabolite, if so split into base name and met #
        met = None
        mat2 = pat2.match(name)
        if mat2 is not None:
            name2 = mat2.group(1)
            met = int(mat2.group(2))
        else: 
            name2 = name
        # locate the compound data structure (CDS)
        cmpd_pickle = find_cmpd_pickle(name2, plate_date)
        if cmpd_pickle is None:
            cmpd_pickle = find_cmpd_pickle(name2, plate_date, sub_commas=False)
        # load the CDS
        with open(cmpd_pickle, 'rb') as pf:
            cds = pickle.load(pf)
        # assemble the necessary info from the CDS
        if met is None:
            smi = cds.SMILES
            monoiso = cds.monoiso
            for ad in cds.adducts:
                if ad.adduct == adduct:
                    d_file = ad.d_file
                    dt = ad.dt
                    mz = ad.monoiso
        else: 
            for mds in cds.metabolites:
                if mds.id == met:
                    monoiso = mds.monoiso
                    inchi = mds.InChI
                    for ad in mds.adducts:
                        if ad.adduct == adduct:
                            d_file = ad.d_file
                            dt = ad.dt
                            mz = ad.monoiso
        # fetch the MS2 spectrum
        spec_txt = spectrum_as_txt(*MassLynxReader(d_file).get_spectrum(3, dt - 0.25, dt + 0.25), mz + 2.)
        # add the spectrum and any metadata into the database
        if met is None:
            cur2.execute(qry2.format(n=plate_n), (cmpd, well, spec_txt, adduct, monoiso, smi))
        else:
            cur2.execute(qry3.format(n=plate_n), (cmpd, well, spec_txt, adduct, monoiso, inchi))

# commit changes to the database and close
con.commit()
con.close()
