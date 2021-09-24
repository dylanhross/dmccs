"""
    Filters out any metabolite data where the MetFrag Score is below a threshold
"""

from sqlite3 import connect
from numpy import array, savetxt

# define queries
# get pre-filter metfrag scores
qry_scores_n = """
SELECT 
    metfrag_score
FROM
    plate_{n}_ms2
JOIN 
    plate_{n}
    ON plate_{n}_ms2.dmim_id = plate_{n}.dmim_id
JOIN
    plate_{n}_id
    ON plate_{n}.dmim_id = plate_{n}_id.dmim_id
WHERE
    -- only metabolites
    met_n > 0
    -- only filter out 'Manual' identifications
    AND notes = 'Manual' 
    -- metfrag score cutoff of 10 determined by parent rank test
    AND metfrag_score IS NOT NULL
;"""

# get dmim_ids to remove from the DMIM DB
qry_dmim_ids_n = """
SELECT 
    plate_{n}.dmim_id
FROM
    plate_{n}
JOIN 
    plate_{n}_ms2
    ON plate_{n}.dmim_id = plate_{n}_ms2.dmim_id
JOIN
    plate_{n}_id
    ON plate_{n}.dmim_id = plate_{n}_id.dmim_id
WHERE
    -- only metabolites
    met_n > 0
    -- only filter out 'Manual' identifications
    AND notes = 'Manual' 
    -- metfrag score cutoff of 10 determined by parent rank test
    AND (
        metfrag_score IS NULL 
        OR metfrag_score < 100
    )
;"""

# get ann_ids associated with the dmim_ids to remove
qry_dmim_annids_n = """
SELECT 
    ann_id
FROM
    plate_{n}_id
WHERE
    dmim_id = ?
;"""

# get str_ids associated with the ann_ids to remove
qry_dmim_strids_n = """
SELECT 
    str_id
FROM
    plate_{n}_3d
WHERE
    ann_id = ?
;"""

# drop from plate_n and plate_n_ms2 tables using dmim_id
qry_drop_dmim_id = [
    """
    DELETE FROM
        plate_{n}
    WHERE
        dmim_id = ?
    ;""",
    """
    DELETE FROM 
        plate_{n}_ms2
    WHERE
        dmim_id = ?
    ;"""
]

# drop from plate_n_id, plate_n_mqn, and plate_n_3d using ann_id
qry_drop_ann_id = [
    """
    DELETE FROM
        plate_{n}_id
    WHERE
        ann_id = ?
    ;""",
    """DELETE FROM
        plate_{n}_mqn
    WHERE
        ann_id = ?
    ;""",
    """DELETE FROM
        plate_{n}_3d
    WHERE
        ann_id = ?
    ;"""
]

# drop from plate_n_md3d using str_id
qry_drop_str_id = """
DELETE FROM
    plate_{n}_md3d
WHERE
    str_id = ?
;"""


def filter_metab_metfrag(dmim_cursor1, dmim_cursor2, dmim_cursor3):
    """ removes all metabolite data for compounds that have bad MS2 scores from MetFrag """
    # dump the pre-filtering metfrag scores
    scores = []
    for n in range(1, 8):
        for score in dmim_cursor1.execute(qry_scores_n.format(n=n)).fetchall():
            scores.append(score)
    savetxt('DMIM_metfrag_scores_prefilter.txt', array(scores))

    # count all the metabolites to remove
    dmim_rem, ann_rem, str_rem = [], [], []
    # iterate through all plates
    for n in range(1, 8):
        for dmim_id in dmim_cursor1.execute(qry_dmim_ids_n.format(n=n)).fetchall():
            dmim_rem.append(dmim_id)
            for ann_id in dmim_cursor2.execute(qry_dmim_annids_n.format(n=n), dmim_id).fetchall():
                ann_rem.append(ann_id)
                for str_id in dmim_cursor3.execute(qry_dmim_strids_n.format(n=n), ann_id).fetchall():
                    str_rem.append(str_id)
    #print(len(dmim_rem), 'dmim_ids to remove by MS2 score')
    print(len(ann_rem), 'ann_ids to remove by MS2 score')
    print(len(str_rem), 'str_ids to remove by MS2 score')
    # drop the selected dmim_ids and ann_ids
    for n in range(1, 8):
        """
        for dmim_id in dmim_rem:
            for qry in qry_drop_dmim_id:
                dmim_cursor1.execute(qry.format(n=n), dmim_id)"""
        for ann_id in ann_rem:
            for qry in qry_drop_ann_id:
                dmim_cursor1.execute(qry.format(n=n), ann_id)
        for str_id in str_rem:
            dmim_cursor1.execute(qry_drop_str_id.format(n=n), str_id)


def main(version):
    """ main execution """
    dmim_fname = 'DMIM_v{version}.db'.format(version=version)

    # initialize the database connections
    dmim_con = connect(dmim_fname)
    dmim_cur1, dmim_cur2, dmim_cur3 = dmim_con.cursor(), dmim_con.cursor(), dmim_con.cursor()

    # transfer parent MS spectra from replicate DB to DMIM DB
    filter_metab_metfrag(dmim_cur1, dmim_cur2, dmim_cur3)

    # commit changes and close DB connections
    dmim_con.commit()
    dmim_con.close()
