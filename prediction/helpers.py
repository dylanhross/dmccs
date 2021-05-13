"""
    Useful helper functions without a good home
"""

from numpy import array, concatenate

from build.add_mqns import smi_to_mqn
from build.add_md3d import compute_3d_descriptors


def featurize(smis, structures, custom_mqns, custom_md3ds):
    """ computes features the same way as in DmimData.data.DMD with the 'custom' kwarg
        smis is a list of SMILES structures
        structures is a list of 3D structures (xyzmq format, text) """
    # generate the complete set of descriptors first
    mqns = array([smi_to_mqn(smi) for smi in smis])
    md3ds = array([compute_3d_descriptors(structure) for structure in structures])
    X = concatenate([mqns, md3ds], axis=1)
    mq2i = {
        'c': 0, 'f': 1, 'cl': 2, 'br': 3, 
        'i': 4, 's': 5, 'p': 6, 'an': 7,
        'cn': 8, 'ao': 9, 'co': 10, 'hac': 11,
        'hbam': 12, 'hba': 13, 'hbdm': 14,
        'hbd': 15, 'neg': 16, 'pos': 17,
        'asb': 18, 'adb': 19, 'atb': 20, 'csb': 21,
        'cdb': 22, 'ctb': 23, 'rbc': 24,
        'asv': 25, 'adv': 26, 'atv': 27, 'aqv': 28,
        'cdv': 29, 'ctv': 30, 'cqv': 31, 'r3': 32,
        'r4': 33, 'r5': 34, 'r6': 35, 'r7': 36,
        'r8': 37, 'r9': 38, 'rg10': 39, 
        'afr': 40, 'bfr': 41
    }
    md2i = {
        'pmi1': 0, 'pmi2': 1, 'pmi3': 2,  
        'rmd02': 3, 'rmd24': 4, 'rmd46': 5, 'rmd68': 6, 'rmd8p': 7
    }
    if custom_mqns == 'all':
        custom_mqns = list(mq2i.keys())
    if custom_md3ds == 'all':
        custom_md3ds = list(md2i.keys())
    mqn_indices = [mq2i[_] for _ in custom_mqns]
    md3d_indices = [md2i[_] for _ in custom_md3ds]
    # shift md3d indices by 42 since they are added at the end 
    md3d_indices = [i + 42 for i in md3d_indices]
    all_indices = mqn_indices + md3d_indices
    # keep only the desired indices
    X = array([_[all_indices] for _ in X])
    return X




