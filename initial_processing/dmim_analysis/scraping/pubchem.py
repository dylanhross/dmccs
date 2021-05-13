"""
    pubchem.py
    Dylan H. Ross
    2018/09/17

    description:
        Tool for using PubChem's REST API to search for compounds and obtain SMILES structures.
        more info here:
            http://pubchemdocs.ncbi.nlm.nih.gov/pug-rest-tutorial
"""


import os
from sqlite3 import connect
from time import sleep
from requests import Session


def name_fetch_cids(session, search_cache_cur, name):
    """
name_fetch_cids
    description:
        Searches for a PubChem CID using a compound name, returning a list of results.
        Returns None if any errors.
    parameters:
        session (requests.Session) -- requests session
        search_cache_cur (sqlite3 connection cursor) -- cursor object for performing search cache queries
        name (str) -- name of the compound to search
    returns:
        (list(int) or None) -- list of PubChem CID(s) matching the search name or None if unsuccessful
"""
    # check the search cache first
    qry = 'SELECT cid FROM pubchem WHERE name=?'
    res = [_ for _ in search_cache_cur.execute(qry, (name,)).fetchall()]
    if res:
        return list(res[0])

    # construct the request URL per REST API documentation
    url_prolog = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
    url_input = "compound/name/"
    url_command = "/cids/"
    url_output = "TXT"
    url = url_prolog + url_input + name + url_command + url_output
    try:
        resp = session.get(url).text
        # check for a failed search response
        if resp.split()[0] == "Status:":
            raise RuntimeError(resp)
        # split response on whitespace
        cids = [int(_) for _ in resp.split()]

        if cids:
            # add an entry to the search cache (the first CID in the results)
            qry = 'INSERT INTO pubchem VALUES (?,?,?,?)'
            search_cache_cur.execute(qry, (cids[0], name, None, None))

        return cids
    except:
        # print("Failed to retrieve PubChem CID for compound: {}".format(name), end=" ")
        return None


def cid_fetch_smiles(session, search_cache_cur, cid, canonical=True):
    """
cid_fetch_smiles
    description:
        Fetches the SMILES structure from the record with the specified PubChem
        CID, returning it as a string. Either the Canonical SMILES or the Isomeric
        SMILES may be fetched (see canonical kwarg). The Default is Canonical
        since that is more likely to be present, however, the Isomeric is preferred
        when available since the overall goal will be to produce 3D structures.
    parameters:
        session (requests.Session) -- requests session
        search_cache_cur (sqlite3 connection cursor) -- cursor object for performing search cache queries
        cid (int) -- PubChem CID
        [canonical (bool)] -- Whether to fetch the canonical SMILES (True) or the
                                isomeric SMILES (False) [optional, default=True]
    returns:
        (str) -- SMILES structure or empty string if no results
"""
    # check the search cache first
    qry = 'SELECT smiles FROM pubchem WHERE cid=? AND smiles IS NOT NULL'
    res = [_ for _ in search_cache_cur.execute(qry, (cid,))]
    if res:
        return res[0][0]

    # construct the request URL per REST API documentation
    url_prolog = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
    url_input = "compound/cid/"
    # option to search for Isomeric SMILES
    stype = "Canonical"
    if not canonical:
        stype = "Isomeric"
    url_command = "/property/{}SMILES/".format(stype)
    url_output = "TXT"
    url = url_prolog + url_input + str(cid) + url_command + url_output
    try:
        resp = session.get(url).text
        # check for a failed search response
        if resp.split()[0] == "Status:":
            raise RuntimeError(resp)
        smi = resp.strip()
        # update entry in the search cache
        qry = 'UPDATE pubchem SET smiles=? WHERE cid=?'
        search_cache_cur.execute(qry, (smi, cid))

        return smi
    except:
        # print("Failed to retrieve {} SMILES for PubChem CID: {}".format(stype, cid), end=" ")
        return ""


def cid_fetch_monoiso(session, search_cache_cur, cid):
    """
cid_fetch_smiles
    description:
        Fetches the monoisotopic mass from the record with the specified PubChem
        CID, returning it as a float.
    parameters:
        session (requests.Session) -- requests session
        search_cache_cur (sqlite3 connection cursor) -- cursor object for performing search cache queries
        cid (int) -- PubChem CID
    returns:
        (float) -- monoisotopic mass or None if no results
"""
    # check the search cache first
    qry = 'SELECT monoiso FROM pubchem WHERE cid=? AND monoiso IS NOT NULL'
    res = [_ for _ in search_cache_cur.execute(qry, (cid,))]
    if res:
        return res[0][0]

    # construct the request URL per REST API documentation
    url_prolog = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
    url_input = "compound/cid/"
    url_command = "/property/MonoisotopicMass/"
    url_output = "TXT"
    url = url_prolog + url_input + str(cid) + url_command + url_output
    try:
        resp = session.get(url).text
        # check for a failed search response
        if resp.split()[0] == "Status:":
            raise RuntimeError(resp)
        monoiso = float(resp.strip())
        # update entry in the search cache
        qry = 'UPDATE pubchem SET monoiso=? WHERE cid=?'
        search_cache_cur.execute(qry, (monoiso, cid))

        return monoiso 
    except:
        # print("Failed to retrieve monoisotopic mass for PubChem CID: {}".format(cid), end=" ")
        return None


def get_pchem_monoiso(name):
    """
get_pchem_monoiso
    description:
        searches PubChem for a monoisotopic mass of a compound by name
    parameters:
        name (str) -- compound name
    returns:
        (float or None) -- monoisotopic mass of the compound (or None if anything goes wrong)
"""
    # start a requests session
    session = Session()

    # load the search cache
    search_cache_con = connect(os.path.join(os.path.split(__file__)[0], 'search_cache/scraping_cache.db'))
    search_cache_cur = search_cache_con.cursor()

    sleep(0.25)
    cids = name_fetch_cids(session, search_cache_cur, name)
    smi = None
    if cids:
        sleep(0.25)
        smi = cid_fetch_monoiso(session, search_cache_cur, cids[0])

    # save the search cache
    search_cache_con.commit()
    search_cache_con.close()

    return smi


def get_pchem_smi(name):
    """
get_pchem_smi
    description:
        searches PubChem for a SMILES structure of a compound by name
    parameters:
        name (str) -- compound name
    returns:
        (str or None) -- SMILES structure of the compound (or None if anything goes wrong)
"""
    # start a requests session
    session = Session()

    # load the search cache
    search_cache_con = connect(os.path.join(os.path.split(__file__)[0], 'search_cache/scraping_cache.db'))
    search_cache_cur = search_cache_con.cursor()

    sleep(0.25)
    cids = name_fetch_cids(session, search_cache_cur, name)
    smi = None
    if cids:
        sleep(0.25)
        smi = cid_fetch_smiles(session, search_cache_cur, cids[0], canonical=False)
        if not smi:
            sleep(0.25)
            smi = cid_fetch_smiles(session, search_cache_cur, cids[0])

        # save the search cache
        search_cache_con.commit()
        search_cache_con.close()

    return smi
