"""
    Adds all of the annotations from the plate_N_id tables
"""

from sqlite3 import connect
from json import load as jload


# define queries
# get parent data from the DMIM DB
qry_dmim_parent_plate_n = """
SELECT
    dmim_id, name
FROM 
    plate_{n}
WHERE
    met_n = 0
;"""

# insert annotation for plate N
qry_dmim_plate_n_id = """
INSERT INTO plate_{n}_id
    (dmim_id, ann_id, annotation, smi, notes)
VALUES
    (?,?,?,?,?)
;"""

# get metabolite data from the DMIM DB
qry_dmim_metab_plate_n = """
SELECT
    dmim_id, name, well, met_n
FROM 
    plate_{n}
WHERE
    met_n > 0
;"""

# get metabolite annotations from the replicate DB
qry_rep_plate_n_id = """
SELECT 
    name, smi, notes
FROM
    plate_{n}_id
WHERE
    well = ?
    AND met_n = ?
;"""

# find dmim_ids that have duplicate annotations
qry_dmim_dup_dmim_ids = """
SELECT 
    COUNT(*) AS c, dmim_id 
FROM 
    plate_{n}_id 
GROUP BY 
    dmim_id, annotation, smi, notes 
    HAVING 
        c > 1
;"""

# select ann_ids of duplicate annotations
qry_dmim_dup_ann_ids = """
SELECT 
    ann_id 
FROM 
    plate_{n}_id 
WHERE
    dmim_id = ?
;"""

# select ann_ids of duplicate annotations
qry_drop_dup_ann_ids = """
DELETE FROM 
    plate_{n}_id 
WHERE
    ann_id = ?
;"""


def parent_annotations(dmim_cursor1, dmim_cursor2):
    """ annotations can be taken directly from the main table for parent compounds
        and SMILES structures taken from parent_data_rmsd.json """

    ann_id = 1

    # load parent_data_rmsd.json
    with open('parent_data_rmsd.json', 'r') as j:
        pdata = jload(j)

    # iterate through all plates
    for n in range(1, 8):
        for dmim_id, name in dmim_cursor1.execute(qry_dmim_parent_plate_n.format(n=n)).fetchall():
            matched = False
            for cmpd in pdata:
                if name == cmpd['name']:
                    smi = cmpd['SMILES']
                    qdata = (dmim_id, ann_id, name, smi, 'ParentDrug')
                    dmim_cursor2.execute(qry_dmim_plate_n_id.format(n=n), qdata)
                    matched = True
                    ann_id += 1
                    break
            if not matched:
                #print('ERROR! no annotation for: {} ({})'.format(name, dmim_id))
                pass
    # return the next available annotation ID
    return ann_id


def transfer_annotations(dmim_cursor1, dmim_cursor2, rep_cursor, ann_id_start):
    """ transfer annotations from the replicate DB to the DMIM DB """
    # annotation id
    ann_id = ann_id_start
    # count the number of dmim_ids without annotations
    noann_p, noann_m = 0, 0
    # iterate through all plates
    for n in range(1, 8):
        for dmim_id, name, well, met_n in dmim_cursor1.execute(qry_dmim_metab_plate_n.format(n=n)).fetchall():
            matched = False
            qdata1 = (well, met_n)
            for name, smi, notes in rep_cursor.execute(qry_rep_plate_n_id.format(n=n), qdata1).fetchall():
                    qdata2 = (dmim_id, ann_id, name, smi, notes)
                    dmim_cursor2.execute(qry_dmim_plate_n_id.format(n=n), qdata2)
                    ann_id += 1
                    matched = True
            if not matched:
                if met_n > 0:
                    noann_m += 1
                else:
                    noann_p += 1
                #print('ERROR! annotation for: {} {} {}'.format(dmim_id, well, met_n))
    print(noann_p + noann_m, 'have no verified annotations', '({} parents, {} metabolites)'.format(noann_p, noann_m))


def remove_duplicate_annotations(dmim_cursor1, dmim_cursor2):
    """ remove any duplicates in the annotations (dmim_id, annotation, smi, and notes are the same) """
    # iterate through all plates
    for n in range(1, 8):
        rm_anns = []
        for c, dmim_id in dmim_cursor1.execute(qry_dmim_dup_dmim_ids.format(n=n)).fetchall():
            rm_anns += [_[0] for _ in dmim_cursor2.execute(qry_dmim_dup_ann_ids.format(n=n), (dmim_id,))][1:]
        for rm_ann in rm_anns:
            dmim_cursor2.execute(qry_drop_dup_ann_ids.format(n=n), (rm_ann,))
            

def main(version):
    """ main execution """
    dmim_fname = 'DMIM_v{version}.db'.format(version=version)
    rep_fname = 'replicate_data_assigned_MS2.db'

    # initialize the database connections
    dmim_con = connect(dmim_fname)
    dmim_cur1, dmim_cur2 = dmim_con.cursor(), dmim_con.cursor()
    rep_con = connect(rep_fname)
    rep_cur = rep_con.cursor()

    # make annotations for parent compounds
    aid = parent_annotations(dmim_cur1, dmim_cur2)

    # transfer data from replicate DB to DMIM DB
    transfer_annotations(dmim_cur1, dmim_cur2, rep_cur, aid)

    # remove duplicate annotations
    remove_duplicate_annotations(dmim_cur1, dmim_cur2)

    # commit changes and close DB connections
    dmim_con.commit()
    dmim_con.close()
    rep_con.close()


