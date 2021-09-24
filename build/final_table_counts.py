"""
    Prints final counts in all of the tables of the DMIM database 
"""

from sqlite3 import connect


# define queries
# count entries in plate_n tables
qry_dmim_plate_n = """
SELECT 
    COUNT(*)
FROM
    plate_{n}
{where}
;"""

# count entries in plate_n_id tables
qry_dmim_plate_n_id = """
SELECT 
    COUNT(*)
FROM
    plate_{n}_id
JOIN
    plate_{n} 
    ON plate_{n}_id.dmim_id = plate_{n}.dmim_id
{where}
;"""

# count entries in plate_n_ms2 tables
qry_dmim_plate_n_ms2 = """
SELECT 
    COUNT(*)
FROM
    plate_{n}_ms2
{where}
;"""

# count entries in plate_n_mqn tables
qry_dmim_plate_n_mqn = """
SELECT 
    COUNT(*)
FROM
    plate_{n}_mqn
{where}
;"""

# count entries in plate_n_3d tables
qry_dmim_plate_n_3d = """
SELECT 
    COUNT(*)
FROM 
    plate_{n}_3d 
    JOIN 
        plate_{n}_id 
        ON plate_{n}_3d.ann_id = plate_{n}_id.ann_id 
    JOIN plate_{n}
        ON plate_{n}_id.dmim_id = plate_{n}.dmim_id
{where}
;"""


# count entries in plate_n_md3d tables
qry_dmim_plate_n_md3d = """
SELECT 
    COUNT(*)
FROM
    plate_{n}_md3d
{where}
;"""


def plate_n_counts(dmim_cursor):
    """ count entries in the plate_n tables """
    print('plate_N\nplate       N    [M+H]+ [M+Na]+ [M+K]+ [M+H-H2O]+ [M]+ parents metabolites')
    msg = 'plate_{} {:6d} {:6d} {:6d}  {:6d}  {:6d}  {:6d}    {:6d}   {:6d}'
    totals = [0, 0, 0, 0, 0, 0, 0, 0]
    for n in range(1, 8):
        n_all = int(dmim_cursor.execute(qry_dmim_plate_n.format(n=n, where='')).fetchall()[0][0])
        n_h = int(dmim_cursor.execute(qry_dmim_plate_n.format(n=n, where='WHERE adduct = "[M+H]+"')).fetchall()[0][0])
        n_na = int(dmim_cursor.execute(qry_dmim_plate_n.format(n=n, where='WHERE adduct = "[M+Na]+"')).fetchall()[0][0])
        n_k = int(dmim_cursor.execute(qry_dmim_plate_n.format(n=n, where='WHERE adduct = "[M+K]+"')).fetchall()[0][0])
        n_h2o = int(dmim_cursor.execute(qry_dmim_plate_n.format(n=n, where='WHERE adduct = "[M+H-H2O]+"')).fetchall()[0][0])
        n_m = int(dmim_cursor.execute(qry_dmim_plate_n.format(n=n, where='WHERE adduct = "[M]+"')).fetchall()[0][0])
        n_par = int(dmim_cursor.execute(qry_dmim_plate_n.format(n=n, where='WHERE met_n = 0')).fetchall()[0][0])
        n_met = int(dmim_cursor.execute(qry_dmim_plate_n.format(n=n, where='WHERE met_n > 0')).fetchall()[0][0])
        print(msg.format(n, n_all, n_h, n_na, n_k, n_h2o, n_m, n_par, n_met))
        totals[0] += n_all
        totals[1] += n_h
        totals[2] += n_na
        totals[3] += n_k
        totals[4] += n_h2o
        totals[5] += n_m
        totals[6] += n_par
        totals[7] += n_met
    print('  total {:6d} {:6d} {:6d}  {:6d}  {:6d}  {:6d}    {:6d}   {:6d}'.format(*totals))
    print()


def plate_n_id_counts(dmim_cursor):
    """ count the entries in the annotation tables """
    print('plate_N_id\nplate       N    parents metabolites')
    msg = 'plate_{} {:6d}    {:6d}   {:6d}'
    totals = [0, 0, 0]
    for n in range(1, 8):
        n_all = int(dmim_cursor.execute(qry_dmim_plate_n_id.format(n=n, where='')).fetchall()[0][0])
        n_par = int(dmim_cursor.execute(qry_dmim_plate_n_id.format(n=n, where='WHERE met_n = 0')).fetchall()[0][0])
        n_met = int(dmim_cursor.execute(qry_dmim_plate_n_id.format(n=n, where='WHERE met_n > 0')).fetchall()[0][0])
        totals[0] += n_all
        totals[1] += n_par
        totals[2] += n_met
        print(msg.format(n, n_all, n_par, n_met))
    print('  total {:6d}    {:6d}   {:6d}'.format(*totals))
    print()


def plate_n_ms2_counts(dmim_cursor):
    """ count the entries in the annotation tables """
    print('plate_N_ms2\nplate       N')
    msg = 'plate_{} {:6d}'
    total = 0
    for n in range(1, 8):
        n_all = int(dmim_cursor.execute(qry_dmim_plate_n_ms2.format(n=n, where='')).fetchall()[0][0])
        total += n_all
        print(msg.format(n, n_all))
    print('  total {:6d} '.format(total))
    print()


def plate_n_mqn_counts(dmim_cursor):
    """ count the entries in the annotation tables """
    print('plate_N_mqn\nplate       N')
    msg = 'plate_{} {:6d}'
    total = 0
    for n in range(1, 8):
        n_all = int(dmim_cursor.execute(qry_dmim_plate_n_mqn.format(n=n, where='')).fetchall()[0][0])
        total += n_all
        print(msg.format(n, n_all))
    print('  total {:6d} '.format(total))
    print()


def plate_n_3d_counts(dmim_cursor):
    """ count the entries in the annotation tables """
    print('plate_N_3d\nplate       N    parents metabolites [M+H]+ [M+Na]+ [M+K]+')
    msg = 'plate_{} {:6d}    {:6d}   {:6d}    {:6d}   {:6d}  {:6d}'
    totals = [0, 0, 0, 0, 0, 0]
    for n in range(1, 8):
        n_all = int(dmim_cursor.execute(qry_dmim_plate_n_3d.format(n=n, where='')).fetchall()[0][0])
        n_par = int(dmim_cursor.execute(qry_dmim_plate_n_3d.format(n=n, where='WHERE met_n = 0')).fetchall()[0][0])
        n_met = int(dmim_cursor.execute(qry_dmim_plate_n_3d.format(n=n, where='WHERE met_n > 0')).fetchall()[0][0])
        n_h = int(dmim_cursor.execute(qry_dmim_plate_n_3d.format(n=n, where='WHERE adduct = "[M+H]+"')).fetchall()[0][0])
        n_na = int(dmim_cursor.execute(qry_dmim_plate_n_3d.format(n=n, where='WHERE adduct = "[M+Na]+"')).fetchall()[0][0])
        n_k = int(dmim_cursor.execute(qry_dmim_plate_n_3d.format(n=n, where='WHERE adduct = "[M+K]+"')).fetchall()[0][0])
        totals[0] += n_all
        totals[1] += n_par
        totals[2] += n_met
        totals[3] += n_h
        totals[4] += n_na
        totals[5] += n_k
        print(msg.format(n, n_all, n_par, n_met, n_h, n_na, n_k))
    print('  total {:6d}    {:6d}   {:6d}    {:6d}   {:6d}  {:6d}'.format(*totals))
    print()


def plate_n_md3d_counts(dmim_cursor):
    """ count the entries in the annotation tables """
    print('plate_N_md3d\nplate       N')
    msg = 'plate_{} {:6d}'
    total = 0
    for n in range(1, 8):
        n_all = int(dmim_cursor.execute(qry_dmim_plate_n_md3d.format(n=n, where='')).fetchall()[0][0])
        total += n_all
        print(msg.format(n, n_all))
    print('  total {:6d} '.format(total))
    print()


def main(version):
    """ main execution """
    dmim_fname = 'DMIM_v{version}.db'.format(version=version)

    # initialize the database connections
    dmim_con = connect(dmim_fname)
    dmim_cur = dmim_con.cursor()

    print()
    plate_n_counts(dmim_cur)
    plate_n_id_counts(dmim_cur)
    plate_n_ms2_counts(dmim_cur)
    plate_n_mqn_counts(dmim_cur)
    plate_n_3d_counts(dmim_cur)
    plate_n_md3d_counts(dmim_cur)

    # commit changes and close DB connections
    dmim_con.commit()
    dmim_con.close()
