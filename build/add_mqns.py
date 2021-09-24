#!/usr/local/Cellar/python@3.9/3.9.0_4/bin/python3 
"""
    Computes MQNs for all annotations
"""

from sqlite3 import connect
from rdkit import Chem
from rdkit.Chem import Descriptors



# define queries
# get data from the DMIM DB
qry_dmim_plate_n_id = """
SELECT 
    ann_id, smi
FROM
    plate_{n}_id
;"""

# insert MQNs for plate N
qry_dmim_plate_n_mqn = """
INSERT INTO plate_{n}_mqn
    (ann_id, mqns)
VALUES
    (?,?)
;"""


def smi_to_mqn(smi):
    """ returns MQNs for an input SMILES structure """
    return Descriptors.rdMolDescriptors.MQNs_(Chem.MolFromSmiles(smi))


def compute_mqns(dmim_cursor1, dmim_cursor2):
    """ compute MQNs for all of the annotations """
    # count the number of dmim_ids without annotations
    nomqns = 0
    # iterate through all plates
    for n in range(1, 8):
        for ann_id, smi in dmim_cursor1.execute(qry_dmim_plate_n_id.format(n=n)).fetchall():
            success = False
            try:
                mqns = smi_to_mqn(smi)
                qdata = (ann_id, ' '.join([str(_) for _ in mqns]))
                dmim_cursor2.execute(qry_dmim_plate_n_mqn.format(n=n), qdata)
                success = True
            except Exception as e:
                print(e)
                nomqns += 1
    print(nomqns, 'annotations have no MQNs')


def main(version):
    """ main execution """
    dmim_fname = 'DMIM_v{version}.db'.format(version=version)

    # initialize the database connections
    dmim_con = connect(dmim_fname)
    dmim_cur1, dmim_cur2 = dmim_con.cursor(), dmim_con.cursor()

    # transfer parent MS spectra from replicate DB to DMIM DB
    compute_mqns(dmim_cur1, dmim_cur2)

    # commit changes and close DB connections
    dmim_con.commit()
    dmim_con.close()

