#

from sqlite3 import connect
from numpy import array
import re

from metfrag2 import fragmenter_rank


def main():
    con = connect('replicate_data_assigned_MS2.db')
    cur = con.cursor()
    pat = re.compile(r'_met[0-9]+')
    cmpds, ranks = [], []
    qry = """SELECT cmpd, well, spectrum, adduct, monoiso, formula, inchi, inchi_key FROM plate_{n}_ms2;"""
    # iterate through all 7 plates
    for n in range(3, 8):
        for cmpd, well, spectrum, adduct, monoiso, formula, inchi, inchi_key in cur.execute(qry.format(n=n)).fetchall():
            if adduct not in ['[M+H]+', '[M+Na]+', '[M+K]+'] or pat.search(cmpd) is not None:
                continue
            try:
                print(cmpd, '(p{}{})'.format(n, well), end=' -> ')
                md = {
                    'adduct': adduct,
                    'formula': formula,
                    'monoiso': monoiso,
                    'inchi': inchi,
                    'inchi_key': inchi_key
                }
                # make mass spectrum into arrays
                m, i = array([[float(_a) for _a in _b.split(' ')] for _b in spectrum.split('\n')[:-1]]).T
                # compute score
                rank = fragmenter_rank(m, i, 1000, md)
                print(rank)
                cmpds.append(cmpd)
                ranks.append(rank)
            except KeyboardInterrupt:
                exit()
            except Exception as e:
                print('FAILED')
                print('!    ', e)
    con.close()
    with open('rank_results2.txt', 'w') as f:
        for cmpd, rank in zip(cmpds, ranks):
            f.write('{} {}\n'.format(cmpd, rank))


if __name__ == '__main__':
    main()
