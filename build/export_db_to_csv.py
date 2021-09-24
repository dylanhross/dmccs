#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3


from sqlite3 import connect


qry = """
SELECT
    plate_{n}.dmim_id, well, annotation, adduct, smi, met_n, mz, ccs_avg, ccs_rsd, ccs_r1, ccs_r2, ccs_r3
FROM
    plate_{n}
    JOIN plate_{n}_id
        ON  plate_{n}.dmim_id = plate_{n}_id.dmim_id
;"""


def main():
    # initialize DB connection
    con = connect('DMIM_v1.1.db')
    cur = con.cursor()
    with open('DMIM_export.csv', 'w') as f:
        f.write('dmim_id,plate_id,annotation,adduct,smi,met_n,mz,ccs_avg,ccs_rsd,ccs_r1,ccs_r2,ccs_r3\n')
        # iterate through all 7 plates and get the data
        for n in range(1, 8):
            for dmim_id, well, annotation, adduct, smi, met_n, mz, ccs_avg, ccs_rsd, ccs_r1, ccs_r2, ccs_r3 in cur.execute(qry.format(n=n)):
                plate_id = 'p{}{}'.format(n, well)
                s = '"{}","{}","{}","{}","{}",{},{:.4f},{:.2f},{:.2f},{:.2f},{:.2f},{:.2f}\n'
                f.write(s.format(dmim_id, plate_id, annotation, adduct, smi, met_n, mz, ccs_avg, ccs_rsd, ccs_r1, ccs_r2, ccs_r3))
    # close database connection
    con.close()


if __name__ == '__main__':
    main()

