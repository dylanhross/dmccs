"""
    dmim_analysis/scraping/cactus.py
    Dylan H. Ross
    2019/10/30

    description:
        tools for using cactus to get images of chemical structures
"""


import os
from sqlite3 import connect
from datetime import datetime
from random import randint
from time import sleep
from requests import Session


def timestamp():
    """
timestamp
    description:
        generates an integer timestamp (14 digits) with the following format:
            YYMMDDHHMMSS{random number 00 to 99}
    returns:
        (int) -- timestamp
"""
    return int(datetime.now().strftime('%y%m%d%H%M%S') + '00') + randint(0, 99)


def structure_to_img(session, search_cache_cur, structure):
    """
structure_to_img
    description:
        Searches cactus using a SMILES structure and returns a png image
        Returns None if any errors.
    parameters:
        session (requests.Session) -- requests session
        search_cache_cur (sqlite3 connection cursor) -- cursor object for performing search cache queries
        structure (str) -- chemical structure (either InChI or SMILES)
    returns:
        (bytes or None) -- png image (bytes) or None if unsuccessful
"""
    # check the search cache first
    qry = 'SELECT img FROM cactus WHERE structure=?'
    res = [_ for _ in search_cache_cur.execute(qry, (structure,))]
    if res:
        return res[0][0]

    # construct the request URL 
    url = 'https://cactus.nci.nih.gov/chemical/structure/{}/image?format=png'.format(structure)
    try:
        sleep(0.25)
        resp = session.get(url)
        if '404' in resp.text:
            return None
        else:
            png = resp.content

            # update entry in the search cache
            tstamp = timestamp()
            qry = 'INSERT INTO cactus VALUES (?,?,?)'
            search_cache_cur.execute(qry, (tstamp, structure, png))

            return png
    except Exception as e:
        #print(e)
        return None


def get_cactus_img(structure):
    """
get_cactus_img
    description:
        searches cactus for an image using a structure
    parameters:
        structure (str) -- chemical structure (in SMILES or InChI format)
    returns:
        (bytes or None) -- png image (bytes) or None if unsuccessful
"""
    # start requests session
    session = Session()

    # load the search cache
    search_cache_con = connect(os.path.join(os.path.split(__file__)[0], 'search_cache/scraping_cache.db'))
    search_cache_cur = search_cache_con.cursor()


    png = None
    png = structure_to_img(session, search_cache_con, structure)

    # save the search cache
    search_cache_con.commit()
    search_cache_con.close()

    return png

