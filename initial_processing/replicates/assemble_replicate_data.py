"""
    Put all of the replicate data into a database
"""


from sqlite3 import connect
from csv import reader
from pickle import load as pload
import os
from glob import glob
from re import search, match
from numpy import mean, std


# initialize connection to the database
con = connect('replicate_data.db')
cur = con.cursor()


# iterate through all 7 plates and add in data
metadata = [
    [1,          2,          3,          4,          5,          6,          7         ], 
    ['20191216', '20191209', '20191209', '20191202', '20191202', '20191115', '20191115'], 
    ['20191217', '20191210', '20191210', '20191203', '20191203', '20191116', '20191116'], 
    ['',         '_p2',      '_p3',      '_p4',      '_p5',      '_p6',      '_p7'     ]
]
for plate_n, date_r2, date_r3, ptag in zip(*metadata): 

    # STEP 1 - create the table for the plate
    qry = """
    CREATE TABLE plate_{} (
        cmpd TEXT NOT NULL,
        well TEXT NOT NULL,
        mz REAL NOT NULL,
        ccs_avg REAL,
        ccs_sd REAL,
        ccs_reps INT NOT NULL, 
        ccs_r1 REAL,
        ccs_r2 REAL,
        ccs_r3 REAL  
    );""".format(plate_n)
    cur.execute(qry)


    # STEP 2 - add all of the replicate 1 data directly from the .csv
    qry = "INSERT INTO plate_{} VALUES (?,?,?,?,?,?,?,?,?);".format(plate_n)
    with open('new_plate_{}.csv'.format(plate_n)) as f:
        rdr = reader(f)
        for cmpd, well, mz, ccs in rdr:
            mz, ccs = float(mz), float(ccs)
            qdata = (cmpd, well, mz, None, None, 1, ccs, None, None)
            cur.execute(qry, qdata)


    # STEP 3 - add replicate 2 data from the pickle files
    qry = "UPDATE plate_{} SET ccs_r2=?, ccs_reps=(ccs_reps + 1) WHERE cmpd=?".format(plate_n)
    for ppath in glob("D:/DMIM_HT/analysis/{}/compounds{}/*.pickle".format(date_r2, ptag)):
        with open(ppath, 'rb') as pf:
            cmpd = search('compounds{}/(.+)[.]pickle'.format(ptag), ppath.replace('\\', '/')).group(1)
            ccs = pload(pf)['ccs']
            cur.execute(qry, (ccs, cmpd))


    # STEP 4 - add replicate 3 data from the pickle files
    qry = "UPDATE plate_{} SET ccs_r3=?, ccs_reps=(ccs_reps + 1) WHERE cmpd=?".format(plate_n)
    for ppath in glob("D:/DMIM_HT/analysis/{}/compounds{}/*.pickle".format(date_r3, ptag)):
        with open(ppath, 'rb') as pf:
            cmpd = search('compounds{}/(.+)[.]pickle'.format(ptag), ppath.replace('\\', '/')).group(1)
            ccs = pload(pf)['ccs']
            cur.execute(qry, (ccs, cmpd))

    # STEP 5 - compute average and SD for compounds with >1 replicate
    qry = "SELECT cmpd, ccs_r1, ccs_r2, ccs_r3 FROM plate_{} WHERE ccs_reps > 1".format(plate_n)
    for cmpd, *ccss in cur.execute(qry).fetchall():
        ccs_all = [_ for _ in ccss if _ is not None]
        m, s = mean(ccs_all), std(ccs_all)
        qry = "UPDATE plate_{} SET ccs_avg=?, ccs_sd=? WHERE cmpd=?".format(plate_n)
        cur.execute(qry, (m, s, cmpd))


# save changes to the database
con.commit()
con.close()

