"""    
    DmimData/db_interface.py
    Dylan H. Ross
    2021/01/15
    
    description:
        Functions for getting data from the DMIM_v?.?.db SQLite3 database 

"""


from sqlite3 import connect
from numpy import array


def qry_mz(db_path):
    """
qry_mz
    description:
        returns the set of measured m/z and CCS values and their corresponding name, adduct, and met_n
    parameters:
        db_path (str) -- path to DMIM_v?.?.db database file
    returns:
        mz, ccs, name, adduct, met_n (numpy.ndarray(float, float, str, str, int))
"""
    con = connect(db_path)
    cur = con.cursor()
    # define the query
    qry = """
        SELECT
            mz, ccs_avg, name, adduct, met_n
        FROM
            plate_{n}
    ;"""
    # accumulate values
    mz, ccs, name, adduct, met_n = [], [], [], [], []
    # iterate through all of the plates
    for n in range(1, 8):
        # fetch data from the database
        for mz_, ccs_, name_, adduct_, met_n_ in cur.execute(qry.format(n=n)):
            mz.append(float(mz_))
            ccs.append(float(ccs_))
            name.append(name_)
            adduct.append(adduct_)
            met_n.append(int(met_n_))
    # close the database and return the data as np.ndarrays
    con.close()
    return array(mz), array(ccs), array(name), array(adduct), array(met_n), array(mz)


def qry_mqn(db_path):
    """
qry_mqn
    description:
        returns the set of MQNs and CCS values and corresponding name/annotation, adduct, and met_n
    parameters:
        db_path (str) -- path to DMIM_v?.?.db database file
    returns:
        mqn, ccs, annotation, adduct, met_n (numpy.ndarray((int,), float, str, str, int))
"""
    con = connect(db_path)
    cur = con.cursor()
    # define the query
    qry = """
        SELECT
            mqns, ccs_avg, annotation, adduct, met_n, mz
        FROM
            plate_{n}_mqn
            JOIN plate_{n}_id
                ON plate_{n}_mqn.ann_id = plate_{n}_id.ann_id
            JOIN plate_{n}
                ON plate_{n}_id.dmim_id = plate_{n}.dmim_id
    ;"""
    # accumulate values
    mqn, ccs, annotation, adduct, met_n, mz = [], [], [], [], [], []
    # iterate through all of the plates
    for n in range(1, 8):
        # fetch data from the database
        for mqn_, ccs_, annotation_, adduct_, met_n_, mz_ in cur.execute(qry.format(n=n)):
            mqn.append([int(_) for _ in mqn_.split()])
            ccs.append(float(ccs_))
            annotation.append(annotation_)
            adduct.append(adduct_)
            met_n.append(int(met_n_))
            mz.append(float(mz_))
    # close the database and return the data as np.ndarrays
    con.close()
    return array(mqn), array(ccs), array(annotation), array(adduct), array(met_n), array(mz)


def qry_md3d(db_path):
    """
qry_md3d
    description:
        returns the set of 3D MDs and CCS values and corresponding name/annotation, adduct, and met_n
    parameters:
        db_path (str) -- path to DMIM_v?.?.db database file
    returns:
        md3d, ccs, annotation, adduct, met_n (numpy.ndarray((float,), float, str, str, int))
"""
    con = connect(db_path)
    cur = con.cursor()
    # define the query
    qry = """
        SELECT
            pmi1, pmi2, pmi3, rmd02, rmd24, rmd46, rmd68, rmd8p, ccs_avg, annotation, adduct, met_n, mz
        FROM
            plate_{n}_md3d
            JOIN plate_{n}_3d
                ON plate_{n}_md3d.str_id = plate_{n}_3d.str_id
            JOIN plate_{n}_id
                ON plate_{n}_3d.ann_id = plate_{n}_id.ann_id
            JOIN plate_{n}
                ON plate_{n}_id.dmim_id = plate_{n}.dmim_id
    ;"""
    # accumulate values
    md3d, ccs, annotation, adduct, met_n, mz = [], [], [], [], [], []
    # iterate through all of the plates
    for n in range(1, 8):
        # fetch data from the database
        for pmi1, pmi2, pmi3, rmd02, rmd24, rmd46, rmd68, rmd8p, ccs_, annotation_, adduct_, met_n_, mz_ in cur.execute(qry.format(n=n)):
            md3d.append([pmi1, pmi2, pmi3, rmd02, rmd24, rmd46, rmd68, rmd8p])
            ccs.append(float(ccs_))
            annotation.append(annotation_)
            adduct.append(adduct_)
            met_n.append(int(met_n_))
            mz.append(float(mz_))
    # close the database and return the data as np.ndarrays
    con.close()
    return array(md3d), array(ccs), array(annotation), array(adduct), array(met_n), array(mz)


def qry_combined(db_path):
    """
qry_combined
    description:
        returns the set of MQNs and MD3Ds and CCS values and corresponding name/annotation, adduct, and met_n
    parameters:
        db_path (str) -- path to DMIM_v?.?.db database file
    returns:
        X, ccs, annotation, adduct, met_n (numpy.ndarray((float,), float, str, str, int))
"""
    con = connect(db_path)
    cur = con.cursor()
    # define the query
    qry = """
        SELECT
            mqns, pmi1, pmi2, pmi3, rmd02, rmd24, rmd46, rmd68, rmd8p, ccs_avg, annotation, adduct, met_n, mz
        FROM
            plate_{n}_md3d
            JOIN plate_{n}_3d
                ON plate_{n}_md3d.str_id = plate_{n}_3d.str_id
            JOIN plate_{n}_id
                ON plate_{n}_3d.ann_id = plate_{n}_id.ann_id
            JOIN plate_{n}_mqn
                ON plate_{n}_id.ann_id = plate_{n}_mqn.ann_id
            JOIN plate_{n}
                ON plate_{n}_id.dmim_id = plate_{n}.dmim_id
    ;"""
    # accumulate values
    mqn, md3d, ccs, annotation, adduct, met_n, mz = [], [], [], [], [], [], []
    # iterate through all of the plates
    for n in range(1, 8):
        # fetch data from the database
        for mqn_, pmi1, pmi2, pmi3, rmd02, rmd24, rmd46, rmd68, rmd8p, ccs_, annotation_, adduct_, met_n_, mz_ in cur.execute(qry.format(n=n)):
            mqn.append([int(_) for _ in mqn_.split()])
            md3d.append([pmi1, pmi2, pmi3, rmd02, rmd24, rmd46, rmd68, rmd8p])
            ccs.append(float(ccs_))
            annotation.append(annotation_)
            adduct.append(adduct_)
            met_n.append(int(met_n_))
            mz.append(float(mz_))
    # combine mqn and md3d into single vector
    x = [a_ + b_ for a_, b_ in zip(mqn, md3d)]
    # close the database and return the data as np.ndarrays
    con.close()
    return array(x), array(ccs), array(annotation), array(adduct), array(met_n), array(mz)

