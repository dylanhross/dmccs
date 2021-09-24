#!/usr/local/Cellar/python@3.9/3.9.0_4/bin/python3 
"""
    Computes MQNs for all annotations
"""

from sqlite3 import connect
from json import load as jload


# define queries
# get data from the DMIM DB
qry_dmim_plate_n_id = """
SELECT 
    ann_id, smi, annotation, adduct
FROM
    plate_{n}_id
    JOIN
        plate_{n}
        ON plate_{n}_id.dmim_id = plate_{n}.dmim_id
;"""

# insert MQNs for plate N
qry_dmim_plate_n_3d = """
INSERT INTO plate_{n}_3d
    (ann_id, str_id, structure)
VALUES
    (?,?,?)
;"""


def get_structures(dmim_cursor1, dmim_cursor2):
    """ annotations can be taken directly from the main table for parent compounds
        and SMILES structures taken from parent_data_rmsd.json """

    # load parent_data_rmsd.json
    with open('parent_data_rmsd.json', 'r') as j:
        pdata = jload(j)
    # load metab_data_fixed_rmsd.json
    with open('metab_data_fixed_rmsd.json', 'r') as j:
        mdata = jload(j)

    # count annotations without structures
    nostr = 0

    # unique structure ID 
    str_id = 1

    # iterate through all plates
    for n in range(1, 8):
        for ann_id, smi, annotation, adduct in dmim_cursor1.execute(qry_dmim_plate_n_id.format(n=n)).fetchall():
            matched = False
            for cmpd in pdata + mdata:
                if annotation == cmpd['name'] and smi == cmpd['SMILES'] and adduct == cmpd['adduct'] and 'structures' in cmpd and cmpd['structures'] != []:
                    matched = True
                    for structure in cmpd['structures']:
                        qdata = (ann_id, str_id, structure)
                        dmim_cursor2.execute(qry_dmim_plate_n_3d.format(n=n), qdata)
                        str_id += 1
                    break
            if not matched:
                nostr += 1
    print(nostr, 'annotations do not have 3D structure(s)')


def main(version):
    """ main execution """
    dmim_fname = 'DMIM_v{version}.db'.format(version=version)

    # initialize the database connections
    dmim_con = connect(dmim_fname)
    dmim_cur1, dmim_cur2 = dmim_con.cursor(), dmim_con.cursor()

    # fetch 3D structures for parent compounds
    get_structures(dmim_cur1, dmim_cur2)

    # commit changes and close DB connections
    dmim_con.commit()
    dmim_con.close()

