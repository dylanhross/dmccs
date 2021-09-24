"""
    Initializes a new database with empty tables
"""

import os
from sqlite3 import connect


# define the table schemas: plate_n, plate_n_id, plate_n_ms2, plate_n_mqn, plate_n_3d
# main tables, contains measurement data
plate_n_schema = """
CREATE TABLE plate_{n} (
    -- global unique identifier
    dmim_id TEXT UNIQUE NOT NULL,
    -- well location on plate
    well TEXT NOT NULL,
    -- analyte name
    name TEXT NOT NULL, 
    -- metabolite number (0 for parent compounds)
    met_n INT NOT NULL,
    -- MS adduct
    adduct TEXT NOT NULL,
    -- m/z
    mz REAL NOT NULL,
    -- CCS (average, RSD %, and individual reps 1-3)
    ccs_avg REAL NOT NULL,
    ccs_rsd REAL NOT NULL,
    ccs_r1 REAL NOT NULL,
    ccs_r2 REAL NOT NULL,
    ccs_r3 REAL NOT NULL
)
;"""

# identification tables, contains metadata corresponding to annotations
plate_n_id_schema = """
CREATE TABLE plate_{n}_id (
    -- global identifier, can have multiple annotations for a single ID (not unique)
    dmim_id TEXT NOT NULL,
    -- annotation identifier, unique
    ann_id INT UNIQUE NOT NULL,
    -- annotation
    annotation TEXT NOT NULL, 
    -- SMILES structure
    smi TEXT NOT NULL,
    -- optional notes
    notes TEXT
)
;"""

# MS2 tables, fragmentation spectra and metfrag score
plate_n_ms2_schema = """
CREATE TABLE plate_{n}_ms2 (
    -- global identifier, only one spectrum for a single ID (unique)
    dmim_id TEXT UNIQUE NOT NULL,
    -- MS2 spectrum
    spectrum TEXT NOT NULL, 
    -- metfrag score, optional only some spectra have these
    metfrag_score REAL
)
;"""

# MQN table, MQNs computed for all annotations
plate_n_mqn_schema = """
CREATE TABLE plate_{n}_mqn (
    -- annotation identifier, only one set of MQNs per annotation (unique)
    ann_id INT UNIQUE NOT NULL,
    -- MQNs
    mqns TEXT NOT NULL
)
;"""

# 3d structure table
plate_n_3d_schema = """
CREATE TABLE plate_{n}_3d (
    -- annotation identifier, often more than 1 structure per annotation (not unique)
    ann_id INT NOT NULL,
    -- structure identifier (unique)
    str_id INT UNIQUE NOT NULL,
    -- 3D structure, xyzmq format (text)
    structure TEXT NOT NULL
)
;"""

# 3d molecular descriptor table
plate_n_md3d_schema = """
CREATE TABLE plate_{n}_md3d (
    -- structure identifier, only 1 set of descriptors per annotation (unique)
    str_id INT UNIQUE NOT NULL,
    -- principal moments of inertia
    pmi1 REAL NOT NULL,
    pmi2 REAL NOT NULL,
    pmi3 REAL NOT NULL,
    -- radial mass distribution
    rmd02 REAL NOT NULL,
    rmd24 REAL NOT NULL,
    rmd46 REAL NOT NULL,
    rmd68 REAL NOT NULL,
    rmd8p REAL NOT NULL
)
;"""


def create_tables(cursor):
    """ creates all of the database tables """
    # iterate through all plates
    for n in range(1, 8):
        cursor.execute(plate_n_schema.format(n=n))
        cursor.execute(plate_n_id_schema.format(n=n))
        cursor.execute(plate_n_ms2_schema.format(n=n))
        cursor.execute(plate_n_mqn_schema.format(n=n))
        cursor.execute(plate_n_3d_schema.format(n=n))
        cursor.execute(plate_n_md3d_schema.format(n=n))


def main(version):
    """ main execution """
    fname = 'DMIM_v{version}.db'.format(version=version)

    # check if the database exists already and remove it if it does
    if os.path.isfile(fname):
        os.remove(fname)

    # initialize the database connection
    con = connect(fname)
    cur = con.cursor()

    # create the tables
    create_tables(cur)

    # commit changes and close DB connection
    con.commit()
    con.close()


