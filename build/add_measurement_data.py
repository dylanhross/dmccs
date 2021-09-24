"""
    Adds all of the CCS measurement data (from replicate data DB)
    Only adds entries which have 3 replicate CCS values and abs(RSD) < 5%
"""

from sqlite3 import connect
import re


# define queries
# fetch data from plate N, only compounds that have 3 CCS reps
qry_rep_plate_n = """
SELECT 
    cmpd, well, mz, ccs_avg, ccs_sd, ccs_r1, ccs_r2, ccs_r3
FROM
    plate_{n}
WHERE
    ccs_reps = 3
;"""

# add an entry to plate N of the DMIM DB
qry_dmim_plate_n = """
INSERT INTO plate_{n}
    (dmim_id, well, name, met_n, adduct, mz, ccs_avg, ccs_rsd, ccs_r1, ccs_r2, ccs_r3)
VALUES
    (?,?,?,?,?,?,?,?,?,?,?)
;"""

# define regex patterns
met_n_pat = re.compile(r'met([0-9]{3})')


def parse_cmpd(cmpd):
    """ takes the compound identifier from the replicate DB and parses out name, adduct, and met_n (if present) 
        returns: name, adduct, met_n """
    cmpd_split = cmpd.split('_')
    adduct = cmpd_split[-1]
    name = '_'.join(cmpd_split[:-1])
    mat = met_n_pat.search(name)
    met_n = int(mat.group(1)) if mat is not None else 0
    return name, met_n, adduct


def sd_to_rsd(ccs_avg, ccs_sd):
    """ converts CCS SD into RSD % """
    return 100. * ccs_sd / ccs_avg 


def transfer_data(dmim_cursor, rep_cursor):
    """ transfer data from replicate DB to DMIM DB """
    # DMIM ID starts at 1
    i = 1
    n_p, n_m = 0, 0

    # iterate through all plates
    for n in range(1, 8):
        # iterate through all of the compounds in the replicate DB
        for cmpd, well, mz, ccs_avg, ccs_sd, ccs_r1, ccs_r2, ccs_r3 in rep_cursor.execute(qry_rep_plate_n.format(n=n)).fetchall():
            name, met_n, adduct = parse_cmpd(cmpd)
            ccs_rsd = sd_to_rsd(ccs_avg, ccs_sd)
            if abs(ccs_rsd) < 5.:
                dmim_id = 'DMIM{:05d}'.format(i)
                qdata = (dmim_id, well, name, met_n, adduct, mz, ccs_avg, ccs_rsd, ccs_r1, ccs_r2, ccs_r3)
                dmim_cursor.execute(qry_dmim_plate_n.format(n=n), qdata)
                i += 1
                if met_n > 0:
                    n_m += 1
                else:
                    n_p += 1
    print('{} entries total'.format(i - 1), '({} parents, {} metabolites)'.format(n_p, n_m))


def main(version):
    """ main execution """
    dmim_fname = 'DMIM_v{version}.db'.format(version=version)
    rep_fname = 'replicate_data_assigned_MS2.db'

    # initialize the database connections
    dmim_con = connect(dmim_fname)
    dmim_cur = dmim_con.cursor()
    rep_con = connect(rep_fname)
    rep_cur = rep_con.cursor()

    # transfer data from replicate DB to DMIM DB
    transfer_data(dmim_cur, rep_cur)

    # commit changes and close DB connections
    dmim_con.commit()
    dmim_con.close()
    rep_con.close()


