#!/usr/local/Cellar/python@3.9/3.9.0_4/bin/python3 
"""
    Computes MQNs for all annotations
"""

from sqlite3 import connect
from io import StringIO
from numpy import loadtxt, array, sum, abs, sqrt, dot, histogram, concatenate
from numpy.linalg import eigh


# define queries
# get data from the DMIM DB
qry_dmim_plate_n_3d = """
SELECT 
    str_id, structure
FROM
    plate_{n}_3d
;"""

# insert MQNs for plate N
qry_dmim_plate_n_md3d = """
INSERT INTO plate_{n}_md3d
    (str_id, pmi1, pmi2, pmi3, rmd02, rmd24, rmd46, rmd68, rmd8p)
VALUES
    (?,?,?,?,?,?,?,?,?)
;"""


def load_and_center_xyzmq(str_data):
    *xyz, m, q = loadtxt(StringIO(str_data), unpack=True)
    xyz = array(xyz)
    tm = sum(m)
    com = array([sum(m * _ / tm) for _ in xyz])
    mcxyz = (xyz.T - com).T
    return mcxyz, m


def inertia_tensor(mcxyz, m):
    """
        mcxyz - centered x, y, and z coords, shape: (3, n_atoms)
        m - atom masses, shape: (n_atoms,)
    """
    cx, cy, cz = mcxyz
    ixx = sum(m * (cy**2. + cz**2.))
    iyy = sum(m * (cx**2. + cz**2.))
    izz = sum(m * (cy**2. + cx**2.))
    ixy = -sum(m * cx * cy)
    ixz = -sum(m * cx * cz)
    iyz = -sum(m * cz * cy)
    I = array(
        [[ixx, ixy, ixz],
         [ixy, iyy, iyz],
         [ixz, iyz, izz]]
    )
    return I


def align_to_principal_axes(mcxyz, m):
    """
        returns rotated coordinates with x, y, and z axes aligned with principal axes
    """
    i_cxyz = inertia_tensor(mcxyz, m)
    eig, axes = eigh(i_cxyz)
    ev0, ev1, ev2 = axes.T
    rxyz = dot(mcxyz.T, axes).T
    return rxyz


def compute_pmi(mcxyz, m):
    rxyz = align_to_principal_axes(mcxyz, m)
    i_mcxyz = inertia_tensor(rxyz, m)
    eig, axes = eigh(i_mcxyz)
    I1, I2, I3 = eig
    return I1, I2, I3


def mass_dist(mcxyz):
    # collect distances of each mass from COM and each charge from COQ
    mdist = []
    for _xyz in mcxyz.T:
        xm, ym, zm = _xyz
        mdist.append(sqrt(xm**2. + ym**2. + zm**2.))
    return mdist


def mass_hist(mdist):
    bins = [0, 2, 4, 6, 8, 20]  # hard-coded based on distributions of all structures
    mhist, _ = histogram(mdist, bins=bins, density=True)
    return mhist


def compute_3d_descriptors(structure):
    # load and center
    mcxyz, m = load_and_center_xyzmq(structure)
    # compute PMI
    I1, I2, I3 = compute_pmi(mcxyz, m)
    # mass and charge distance distributions, binned
    mhist = mass_hist(mass_dist(mcxyz))
    # put everything together and return as an array
    return I1, I2, I3, *mhist


def add_md3d(dmim_cursor1, dmim_cursor2):
    """ compute 3D molecular descriptors for all 3D structures """
    # count the number of dmim_ids without annotations
    nomd3ds = 0
    # iterate through all plates
    for n in range(1, 8):
        for str_id, structure in dmim_cursor1.execute(qry_dmim_plate_n_3d.format(n=n)).fetchall():
            success = False
            try:
                qdata = (str_id, *compute_3d_descriptors(structure))
                dmim_cursor2.execute(qry_dmim_plate_n_md3d.format(n=n), qdata)
                success = True
            except Exception as e:
                print(e)
                nomd3ds += 1
    print(nomd3ds, 'structures have no MD3Ds')


def main(version):
    """ main execution """
    dmim_fname = 'DMIM_v{version}.db'.format(version=version)

    # initialize the database connections
    dmim_con = connect(dmim_fname)
    dmim_cur1, dmim_cur2 = dmim_con.cursor(), dmim_con.cursor()

    # transfer parent MS spectra from replicate DB to DMIM DB
    add_md3d(dmim_cur1, dmim_cur2)

    # commit changes and close DB connections
    dmim_con.commit()
    dmim_con.close()

