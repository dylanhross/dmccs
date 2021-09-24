"""
    Adds MS2 spectra from the plate_N_ms2 tables
"""

from sqlite3 import connect


# define queries
# get data from the DMIM DB
qry_dmim_plate_n = """
SELECT 
    dmim_id, well, name, adduct
FROM
    plate_{n}
;"""

# get MS2 data from the replicate DB
qry_rep_plate_n_ms2 = """
SELECT 
    spectrum, metfrag_score
FROM
    plate_{n}_ms2
WHERE 
    cmpd = ?
    AND well = ?
    AND adduct = ?
;"""

# insert MS2 for plate N
qry_dmim_plate_n_ms2 = """
INSERT INTO plate_{n}_ms2
    (dmim_id, spectrum, metfrag_score)
VALUES
    (?,?,?)
;"""


def transfer_spectra(dmim_cursor1, dmim_cursor2, rep_cursor):
    """ transfer parent MS2 spectra replicate DB to DMIM DB """
    # count the number of dmim_ids without annotations
    noms2 = 0
    # count null metfrag scores
    score_null_p, score_null_m = 0, 0
    # iterate through all plates
    for n in range(1, 8):
        for dmim_id, well, name, adduct in dmim_cursor1.execute(qry_dmim_plate_n.format(n=n)).fetchall():
            matched = False
            cmpd = '{}_{}'.format(name, adduct)
            qdata1 = (cmpd, well, adduct)
            for spectrum, metfrag_score in rep_cursor.execute(qry_rep_plate_n_ms2.format(n=n), qdata1).fetchall():
                    qdata2 = (dmim_id, spectrum, metfrag_score)
                    dmim_cursor2.execute(qry_dmim_plate_n_ms2.format(n=n), qdata2)
                    matched = True
                    if metfrag_score is None:
                        if '_met' in name:
                            score_null_p += 1
                        else:
                            score_null_m += 1
            if not matched:
                noms2 += 1
                #print('ERROR! annotation for: {} {} {}'.format(dmim_id, well, met_n))
    print(noms2, 'have no MS2 spectra')
    print(score_null_p + score_null_m, 'have NULL metfrag_score', '({} parents, {} metabolites)'.format(score_null_p, score_null_m))


def main(version):
    """ main execution """
    dmim_fname = 'DMIM_v{version}.db'.format(version=version)
    rep_fname = 'replicate_data_assigned_MS2.db'

    # initialize the database connections
    dmim_con = connect(dmim_fname)
    dmim_cur1, dmim_cur2 = dmim_con.cursor(), dmim_con.cursor()
    rep_con = connect(rep_fname)
    rep_cur = rep_con.cursor()

    # transfer parent MS spectra from replicate DB to DMIM DB
    transfer_spectra(dmim_cur1, dmim_cur2, rep_cur)

    # commit changes and close DB connections
    dmim_con.commit()
    dmim_con.close()
    rep_con.close()


