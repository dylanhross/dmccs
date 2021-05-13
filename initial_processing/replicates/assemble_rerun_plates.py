"""
    Put all of the replicate data into a database
"""


from sqlite3 import connect


# initialize connection to the database
con = connect('replicate_data.db')
cur = con.cursor()


# iterate through all 7 plates and aggregate the entries that do not have all replicates
metadata = [
    [1,          2,          3,          4,          5,          6,          7         ], 
    ['20191216', '20191209', '20191209', '20191202', '20191202', '20191115', '20191115'], 
    ['20191217', '20191210', '20191210', '20191203', '20191203', '20191116', '20191116'], 
    ['',         '_p2',      '_p3',      '_p4',      '_p5',      '_p6',      '_p7'     ]
]
for plate_n, date_r2, date_r3, ptag in zip(*metadata): 

    # find compounds where rep2 CCS is NULL
    qry = "SELECT cmpd, well, mz, ccs_r1 FROM plate_{} WHERE ccs_r2 IS NULL;".format(plate_n)
    with open('rerun_p{}r2.csv'.format(plate_n), 'w') as out:
        for cmpd, well, mz, ccs in cur.execute(qry).fetchall():
            mz, ccs = float(mz), float(ccs)
            line = '"{}",{},{:.4f},{:.2f}\n'.format(cmpd, well, mz, ccs)
            out.write(line)

    # find compounds where rep2 CCS is NULL
    qry = "SELECT cmpd, well, mz, ccs_r1 FROM plate_{} WHERE ccs_r3 IS NULL;".format(plate_n)
    with open('rerun_p{}r3.csv'.format(plate_n), 'w') as out:
        for cmpd, well, mz, ccs in cur.execute(qry).fetchall():
            mz, ccs = float(mz), float(ccs)
            line = '"{}",{},{:.4f},{:.2f}\n'.format(cmpd, well, mz, ccs)
            out.write(line)


con.close()

